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

def add_ipe(qc, qubit, cbits, element):
    """Append iterative phase estimation of chi(element) to len(cbits) bits.
    Bit b_k (phi = 0.b1 b2 ... bm) goes into cbits[k-1]."""
    m = len(cbits)
    for k in range(m, 0, -1):
        qc.reset(qubit)
        qc.h(qubit)
        qc.append(oracle_gate(element, 2 ** (k - 1)), [qubit])
        for l in range(k + 1, m + 1):
            with qc.if_test((cbits[l - 1], 1)):
                qc.p(-math.pi / 2 ** (l - k), qubit)
        qc.h(qubit)
        qc.measure(qubit, cbits[k - 1])


def ipe_circuit(elements, m):
    """Read chi(e) to m bits for every e in `elements`, all in ONE shot (same chi)."""
    q = QuantumRegister(1, "q")
    c = ClassicalRegister(m * len(elements), "c")
    qc = QuantumCircuit(q, c)
    for j, el in enumerate(elements):
        add_ipe(qc, q[0], [c[j * m + i] for i in range(m)], el)
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

def is_identity(oracle, terms, max_shots=16, batch=4):
    """Is the product of element**power over `terms` equal to e?"""
    qc = identity_test_circuit(terms)
    ones = shots = 0
    while shots < max_shots:
        counts = oracle.run(qc, shots=batch)
        ones += counts.get("1", 0)
        shots += batch
        if ones >= 2:
            return False
    return True


def reduce_to_order(oracle, element, N):
    for p in sorted(factorize(N)):
        while N % p == 0 and is_identity(oracle, [(element, N // p)]):
            N //= p
    return N


def find_order(oracle, element="g", max_shots=8):
    m = 2 * MAX_BITS + 1
    qc = ipe_circuit([element], m)
    candidates, failed = set(), []
    for i in range(max_shots):
        if not can_afford(oracle, qc):
            log("  out of budget for order finding")
            break
        phi = read_phases(oracle.run(qc, shots=1), [element], m)[element]
        r = phi.limit_denominator(2 ** MAX_BITS).denominator
        new = {r} | {lcm(C, r) for C in candidates}
        new = {N for N in new if N <= 2 ** MAX_BITS} - candidates
        candidates |= new
        log(f"  shot {i}: r = {r}, {len(new)} new candidate(s)")
        for N in sorted(new, reverse=True):
            if any(F % N == 0 for F in failed):
                continue
            if is_identity(oracle, [(element, N)]):
                order = reduce_to_order(oracle, element, N)
                log(f"  g^{N} = e passed; after trimming, order = {order}")
                return order
            failed.append(N)
    return None


def find_exponents(oracle, n, elements=("A", "B", "C"), max_shots=6):
    """Find a with X = g^a for each X, verified on the device.

    Returns (found, forged): found = {X: a}; forged = an element judged forged, or None.
    A 'strike' is a shot whose data contradicts X being a power of g: an off-grid
    reading, an impossible congruence, or a full-modulus candidate that fails verification.
    """
    m = n.bit_length() + MARGIN_BITS
    cands = {X: set() for X in elements}
    failed = {X: set() for X in elements}
    strikes = {X: 0 for X in elements}
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

        x, off = {}, {}
        for name, phi in phases.items():
            y = phi * n
            x[name] = round(y) % n
            off[name] = abs(y - round(y)) > Fraction(1, 4)
        if off["g"]:
            log(f"  joint shot {i}: chi(g) off the 1/n grid -> junk shot, skipped")
            continue

        for X in todo:
            strike = False
            cong = None if off[X] else congruence(x["g"], x[X], n)
            if cong is None:
                strike = True
                log(f"  joint shot {i}: {X} off-grid or impossible congruence")
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
                    return found, X
    return found, None


# ------------------------------------------------------------------ entry

def solve(oracle):
    order = find_order(oracle)
    if order is None:
        log("order finding failed; guessing 1")
        return {"order": 1, "secret": "none"}
    found, forged = find_exponents(oracle, order)
    if forged is None and all(X in found for X in ("A", "B", "C")):
        secret = found["A"] * found["B"] * found["C"] % order
    else:
        if forged is None:
            missing = [X for X in ("A", "B", "C") if X not in found]
            log(f"  never verified: {missing} -> treating as forgery")
        secret = "none"
    log(f"order = {order}, secret = {secret}, queries used = {oracle.queries_used}")
    return {"order": order, "secret": secret}