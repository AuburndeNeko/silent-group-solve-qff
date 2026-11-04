"""Silent Group solver.

Step 3: robust order finding.  Step 4: secret abc mod ord(g) + forgery detection.
"""

import math
import random
from fractions import Fraction

from qiskit import ClassicalRegister, QuantumCircuit, QuantumRegister

from silentgroup import oracle_gate

MAX_BITS = 60       # |G| <= 2^60, so every element's order is <= 2^60
MARGIN_BITS = 4     # extra bits when rounding to the 1/n grid
FORGERY_STRIKES = 3  # contradicting shots before calling an element forged
VERBOSE = True      # set False before submitting


def log(*args):
    if VERBOSE:
        print(*args)

def can_afford(oracle, qc, shots=1):
    """True if this run fits in the remaining budget (None = no limit)."""
    left = oracle.budget_left()
    return left is None or oracle.cost(qc, shots) <= left


# ---------------------------------------------------------------- circuits

def add_ipe(qc, qubit, cbits, element, mult=1):
    """Iterative phase estimation of chi(element**mult) to len(cbits) bits."""
    m = len(cbits)
    for k in range(m, 0, -1):
        qc.reset(qubit)
        qc.h(qubit)
        qc.append(oracle_gate(element, mult * 2 ** (k - 1)), [qubit])
        for l in range(k + 1, m + 1):
            with qc.if_test((cbits[l - 1], 1)):
                qc.p(-math.pi / 2 ** (l - k), qubit)
        qc.h(qubit)
        qc.measure(qubit, cbits[k - 1])


def ipe_circuit(elements, m, mult=1):
    q = QuantumRegister(1, "q")
    c = ClassicalRegister(m * len(elements), "c")
    qc = QuantumCircuit(q, c)
    for j, el in enumerate(elements):
        add_ipe(qc, q[0], [c[j * m + i] for i in range(m)], el, mult)
    return qc


def read_phases(counts, elements, m):
    """Counts from a 1-shot run -> {element: exact Fraction phase}. Never floats."""
    (key,) = counts
    bits = key[::-1]                                  # bits[i] = classical bit i
    return {el: Fraction(int(bits[j * m:(j + 1) * m], 2), 2 ** m)
            for j, el in enumerate(elements)}


def identity_test_circuit(terms):
    """terms = [(element, power), ...]. Always measures 0 iff the product is e."""
    qc = QuantumCircuit(1, 1)
    qc.h(0)
    for element, power in terms:
        qc.append(oracle_gate(element, power), [0])
    qc.h(0)
    qc.measure(0, 0)
    return qc


# ---------------------------------------------------------- number theory

def lcm(a, b):
    return a // math.gcd(a, b) * b


_SMALL_PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]


def is_prime(n):
    if n < 2:
        return False
    for p in _SMALL_PRIMES:
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in _SMALL_PRIMES:
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def _pollard_rho(n):
    while True:
        c = random.randrange(1, n)
        x = y = random.randrange(2, n)
        d = 1
        while d == 1:
            x = (x * x + c) % n
            y = (y * y + c) % n
            y = (y * y + c) % n
            d = math.gcd(abs(x - y), n)
        if d != n:
            return d


def factorize(n):
    out = {}
    for p in range(2, 1000):
        while n % p == 0:
            out[p] = out.get(p, 0) + 1
            n //= p
    stack = [n] if n > 1 else []
    while stack:
        x = stack.pop()
        if is_prime(x):
            out[x] = out.get(x, 0) + 1
        else:
            d = _pollard_rho(x)
            stack += [d, x // d]
    return out


def congruence(xg, xX, n):
    """From xX = a*xg (mod n), return (r, mod) meaning a = r (mod mod), or None if impossible."""
    d = math.gcd(xg, n)
    if xX % d:
        return None
    mod = n // d
    if mod == 1:
        return (0, 1)
    return ((xX // d) * pow(xg // d, -1, mod) % mod, mod)


def crt_merge(c1, c2):
    """Combine a = r1 (mod m1) and a = r2 (mod m2). None if they contradict."""
    (r1, m1), (r2, m2) = c1, c2
    g = math.gcd(m1, m2)
    if (r2 - r1) % g:
        return None
    l = m1 // g * m2
    m2g = m2 // g
    t = 0 if m2g == 1 else ((r2 - r1) // g) * pow(m1 // g, -1, m2g) % m2g
    return ((r1 + m1 * t) % l, l)


# ---------------------------------------------------- quantum subroutines

def is_identity(oracle, terms, max_shots=16, batch=4, reject_at=2):
    """Is the product of element**power over `terms` equal to e?
    Says 'no' once it sees `reject_at` ones (2 normally; 4 when 'no' must be a proof)."""
    qc = identity_test_circuit(terms)
    ones = shots = 0
    while shots < max_shots:
        counts = oracle.run(qc, shots=batch)
        ones += counts.get("1", 0)
        shots += batch
        if ones >= reject_at:
            return False
    return True


def reduce_to_order(oracle, element, N, primes=None):
    """N is a multiple of the order; strip unneeded prime factors (only `primes` if given)."""
    for p in sorted(primes if primes is not None else factorize(N)):
        while N % p == 0 and is_identity(oracle, [(element, N // p)]):
            N //= p
    return N


def find_order(oracle, element="g", max_shots=8, trim=True):
    """Order of `element` (or an untrimmed multiple if trim=False).

    After a candidate N0 fails g^N0 = e, it is often a correct-but-partial divisor of the
    order. Then h = g^N0 has small order k <= 2^MAX_BITS / N0, so reading h's order needs
    far fewer bits than a fresh 121-bit readout. We use that 'leftover' readout whenever
    it is cheaper; if N0 was junk it simply fails and we fall back to a full readout.
    """
    m_full = 2 * MAX_BITS + 1
    qc_full = ipe_circuit([element], m_full)
    candidates, failed = set(), []
    leftover_tries = {}

    def consider(r):
        new = {r} | {lcm(C, r) for C in candidates}
        new = {N for N in new if N <= 2 ** MAX_BITS} - candidates
        candidates.update(new)
        for N in sorted(new, reverse=True):
            if any(F % N == 0 for F in failed):
                continue
            if is_identity(oracle, [(element, N)]):
                return N
            failed.append(N)
        return None

    for i in range(max_shots):
        qc, mode = qc_full, "full"
        if failed:
            N0 = max(failed)
            K = (2 ** MAX_BITS) // N0  # leftover order is at most K
            if K >= 2 and leftover_tries.get(N0, 0) < 2:
                m_h = 2 * K.bit_length() + 2
                qc_h = ipe_circuit([element], m_h, mult=N0)
                if oracle.cost(qc_h, 1) < oracle.cost(qc_full, 1):
                    qc, mode = qc_h, "leftover"
                    leftover_tries[N0] = leftover_tries.get(N0, 0) + 1
        if not can_afford(oracle, qc):
            log("  out of budget for order finding")
            break
        if mode == "full":
            phi = read_phases(oracle.run(qc, shots=1), [element], m_full)[element]
            r = phi.limit_denominator(2 ** MAX_BITS).denominator
            log(f"  121-bit shot {i}: r = {r}")
        else:
            phi = read_phases(oracle.run(qc, shots=1), [element], m_h)[element]
            k = phi.limit_denominator(K).denominator
            r = N0 * k
            log(f"  leftover shot {i} ({m_h} bits on g^{N0}): k = {k} -> try {r}")
            if k == 1:
                leftover_tries[N0] = 2   # g^N0 != e yet its order reads as 1 -> N0 is junk
        N = consider(r)
        if N:
            order = reduce_to_order(oracle, element, N) if trim else N
            log(f"  g^{N} = e passed" + (f"; order = {order}" if trim else " (trimming deferred)"))
            return order
    return None


def find_exponents(oracle, N, elements=("A", "B", "C"), max_shots=6):
    """N: a multiple of ord(g) that passed g^N = e (not yet trimmed).

    Returns (found, forged, n): exponents {X: a}, a forged element or None, and n = ord(g).
    Trimming uses the first clean chi(g) reading: only primes of gcd(x_g, N) can be extra.
    """
    n, trimmed = N, False
    m = N.bit_length() + MARGIN_BITS
    cands = {X: set() for X in elements}
    failed = {X: set() for X in elements}
    strikes = {X: 0 for X in elements}
    power_checked = set()
    found = {}
    for i in range(max_shots):
        todo = [X for X in elements if X not in found]
        if not todo:
            break
        names = ["g"] + todo
        qc = ipe_circuit(names, m)
        if not can_afford(oracle, qc):
            log("  out of budget for exponent finding")
            break
        phases = read_phases(oracle.run(qc, shots=1), names, m)

        if not trimmed:                                   # trim using this shot's chi(g)
            y = phases["g"] * N
            if abs(y - round(y)) > Fraction(1, 4):
                log(f"  joint shot {i}: chi(g) off the 1/N grid -> junk shot, skipped")
                continue
            d = math.gcd(round(y) % N, N)
            if d > 1:
                n = reduce_to_order(oracle, "g", N, primes=factorize(d))
            trimmed = True
            log(f"  order = {n} (gcd(x_g, N) = {d}; only its primes needed checking)")

        x, off = {}, {}
        for name, phi in phases.items():
            y = phi * n
            x[name] = round(y) % n
            off[name] = abs(y - round(y)) > Fraction(1, 4)
        if off["g"]:
            log(f"  joint shot {i}: chi(g) off the 1/n grid -> junk shot, skipped")
            continue

        for X in todo:
            if off[X]:
                if X not in power_checked:
                    power_checked.add(X)
                    if not is_identity(oracle, [(X, n)], reject_at=4):
                        log(f"  {X}^n != e -> {X} is not a power of g: forged (proof)")
                        return found, X, n
                log(f"  joint shot {i}: {X} off-grid but {X}^n = e -> junk reading, skipped")
                continue
            strike = False
            cong = congruence(x["g"], x[X], n)
            if cong is None:
                strike = True
                log(f"  joint shot {i}: {X} gives an impossible congruence")
            else:
                merged = {cong} | {crt_merge(c, cong) for c in cands[X]}
                merged.discard(None)
                cands[X] |= merged
                log(f"  joint shot {i}: {X}: a = {cong[0]} mod {cong[1]}")
                for r in sorted(r for r, mod in merged if mod == n):
                    if r in failed[X]:
                        strike = True
                        continue
                    terms = [(X, 1)] + ([("g", -r)] if r else [])
                    if is_identity(oracle, terms):
                        found[X] = r
                        log(f"  verified: {X} = g^{r}")
                        break
                    failed[X].add(r)
                    strike = True
            if X in found:
                continue
            if strike:
                strikes[X] += 1
                log(f"  {X}: strike {strikes[X]}/{FORGERY_STRIKES}")
                if strikes[X] >= FORGERY_STRIKES:
                    log(f"  {X} contradicts every exponent -> forged; stopping early")
                    return found, X, n
    if not trimmed:                                       # never got a clean chi(g): full trim
        n = reduce_to_order(oracle, "g", N)
    return found, None, n


def solve(oracle):
    N = find_order(oracle, trim=False)
    if N is None:
        log("order finding failed; guessing 1")
        return {"order": 1, "secret": "none"}
    found, forged, order = find_exponents(oracle, N)
    if forged is None and all(X in found for X in ("A", "B", "C")):
        secret = found["A"] * found["B"] * found["C"] % order
    else:
        if forged is None:
            missing = [X for X in ("A", "B", "C") if X not in found]
            log(f"  never verified: {missing} -> treating as forgery")
        secret = "none"
    log(f"order = {order}, secret = {secret}, queries used = {oracle.queries_used}")
    return {"order": order, "secret": secret}