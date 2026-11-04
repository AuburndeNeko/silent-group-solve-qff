"""Silent Group solver v4.1.

Order:     121-bit iterative phase estimation of chi(g), continued fractions, lcm of
           subsets of readings, device test g^N = e. Trimming is deferred: the first clean
           joint-shot reading x_g shows which primes (those of gcd(x_g, N)) could be extra.
Exponents: joint shots read chi(g), chi(A), chi(B), chi(C) in one shot (same character),
           round to the 1/n grid, solve a*x_g = x_X (mod n), merge with CRT, verify on device.
           Partial congruences (a known mod n/d, small d) are finished by testing the d lifts.
Forgery:   proof via X^n != e when X reads off-grid; otherwise 3 strikes (contradicting
           shots). A shot where every element fails at once is blamed on a junk chi(g).
"""

import math
import random
from fractions import Fraction

from qiskit import ClassicalRegister, QuantumCircuit, QuantumRegister

from silentgroup import oracle_gate

MAX_BITS = 60        # |G| <= 2^60, so every element's order is <= 2^60
MARGIN_BITS = 4      # extra bits when rounding to the 1/n grid
FORGERY_STRIKES = 3  # contradicting shots before calling an element forged
COMPLETE_LIMIT = 32  # finish a partial congruence by brute force if n/mod <= this
VERBOSE = True       # set False before submitting
ELEMENTS = ("A", "B", "C")


def log(*args):
    if VERBOSE:
        print(*args)


def can_afford(oracle, qc, shots=1):
    """True if this run fits in the remaining budget (None = no limit)."""
    left = oracle.budget_left()
    return left is None or oracle.cost(qc, shots) <= left


# ---------------------------------------------------------------- circuits

def add_ipe(qc, qubit, cbits, element):
    """Iterative phase estimation of chi(element) to len(cbits) bits.
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
    """Deterministic Miller-Rabin, valid far beyond 2^64."""
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
    """Prime factorization as {prime: exponent}."""
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
    """From xX = a*xg (mod n): (r, mod) meaning a = r (mod mod), or None if impossible."""
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
    Says 'no' once it sees `reject_at` ones (2 normally; 4 when 'no' must be a proof;
    1 for a cheap filter). Out of budget -> 'no' (never accept unchecked)."""
    qc = identity_test_circuit(terms)
    ones = shots = 0
    while shots < max_shots:
        if not can_afford(oracle, qc, batch):
            log("  out of budget during an identity test -> treating as 'no'")
            return False
        counts = oracle.run(qc, shots=batch)
        ones += counts.get("1", 0)
        shots += batch
        if ones >= reject_at:
            return False
    return True


def verify(oracle, X, a):
    """Full test: is X == g^a?"""
    return is_identity(oracle, [(X, 1)] + ([("g", -a)] if a else []))


def complete_partial(oracle, X, r, mod, n, limit=COMPLETE_LIMIT):
    """We know a = r (mod mod) with mod | n. Try each of the n/mod lifts with a cheap
    filter test, then fully verify the survivor. Returns a, or None."""
    d = n // mod
    if d > limit:
        return None
    for j in range(d):
        a = (r + j * mod) % n
        terms = [(X, 1)] + ([("g", -a)] if a else [])
        if is_identity(oracle, terms, max_shots=8, batch=4, reject_at=1) and verify(oracle, X, a):
            return a
    return None


def reduce_to_order(oracle, element, N, primes=None):
    """N is a multiple of the order; strip unneeded prime factors (only `primes` if given)."""
    for p in sorted(primes if primes is not None else factorize(N)):
        while N % p == 0 and is_identity(oracle, [(element, N // p)]):
            N //= p
    return N


def find_order(oracle, element="g", max_shots=8, trim=True):
    """Order of `element` (or an untrimmed multiple of it if trim=False), robust to junk.

    If a failed candidate N0 is backed by >= 2 readings (readings dividing it), it is
    almost surely a correct-but-partial divisor of the order (junk never repeats), so we
    try N0*d for small d before spending another 121-bit readout.
    """
    m = 2 * MAX_BITS + 1
    qc = ipe_circuit([element], m)
    candidates, failed, readings = set(), [], []
    tried_lift = set()

    def done(N):
        order = reduce_to_order(oracle, element, N) if trim else N
        log(f"  g^{N} = e passed" + (f"; order = {order}" if trim else " (trimming deferred)"))
        return order

    for i in range(max_shots):
        if not can_afford(oracle, qc):
            log("  out of budget for order finding")
            break
        phi = read_phases(oracle.run(qc, shots=1), [element], m)[element]
        r = phi.limit_denominator(2 ** MAX_BITS).denominator
        readings.append(r)
        new = {r} | {lcm(C, r) for C in candidates}
        new = {N for N in new if N <= 2 ** MAX_BITS} - candidates
        candidates |= new
        log(f"  121-bit shot {i}: r = {r}, {len(new)} new candidate(s)")
        for N in sorted(new, reverse=True):
            if any(F % N == 0 for F in failed):
                continue
            if is_identity(oracle, [(element, N)]):
                return done(N)
            failed.append(N)

        supported = [N for N in failed if N not in tried_lift
                     and sum(1 for x in readings if x > 1 and N % x == 0) >= 2]
        if supported:
            N0 = max(supported)
            tried_lift.add(N0)
            log(f"  {N0} is backed by >= 2 readings -> trying small multiples")
            for d in range(2, COMPLETE_LIMIT + 1):
                if N0 * d > 2 ** MAX_BITS:
                    break
                if (is_identity(oracle, [(element, N0 * d)], max_shots=8, batch=4, reject_at=1)
                        and is_identity(oracle, [(element, N0 * d)])):
                    log(f"  lifted partial order: {N0} x {d}")
                    return done(N0 * d)
    return None


def find_exponents(oracle, N, elements=ELEMENTS, max_shots=6, hard_cap=12):
    """N: a multiple of ord(g) that passed g^N = e (not yet trimmed).

    Returns (found, forged, n): exponents {X: a}, a forged element or None, n = ord(g).
    After `max_shots` joint shots, only strike-free elements keep going (up to hard_cap).
    """
    n, trimmed = N, False
    m = N.bit_length() + MARGIN_BITS
    cands = {X: set() for X in elements}
    failed = {X: set() for X in elements}
    strikes = {X: 0 for X in elements}
    power_checked, tried_complete = set(), set()
    found = {}

    for i in range(hard_cap):
        todo = [X for X in elements if X not in found]
        if not todo:
            break
        if i >= max_shots and all(strikes[X] > 0 for X in todo):
            break                                         # only suspects left
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

        shot_strikes, outcome = [], {}
        for X in todo:
            if off[X]:
                if X not in power_checked:
                    power_checked.add(X)
                    if not is_identity(oracle, [(X, n)], reject_at=4):
                        log(f"  {X}^n != e -> {X} is not a power of g: forged (proof)")
                        return found, X, n
                log(f"  joint shot {i}: {X} off-grid but {X}^n = e -> junk reading, skipped")
                outcome[X] = "skipped"
                continue
            strike = failed_full = False
            cong = congruence(x["g"], x[X], n)
            if cong is None:
                strike = True
                outcome[X] = "impossible"
                log(f"  joint shot {i}: {X} gives an impossible congruence")
            else:
                merged = {cong} | {crt_merge(c, cong) for c in cands[X]}
                merged.discard(None)
                cands[X] |= merged
                log(f"  joint shot {i}: {X}: a = {cong[0]} mod {cong[1]}")
                for r in sorted(r for r, mod in merged if mod == n):
                    if r in failed[X]:
                        strike = failed_full = True
                        continue
                    if verify(oracle, X, r):
                        found[X] = r
                        log(f"  verified: {X} = g^{r}")
                        break
                    failed[X].add(r)
                    strike = failed_full = True
                outcome[X] = ("verified" if X in found
                              else "failed_full" if failed_full else "partial")
            if X not in found and strike:
                shot_strikes.append(X)

        for X in todo:                                    # finish partial congruences
            if X in found or not cands[X]:
                continue
            r, mod = max(cands[X], key=lambda c: c[1])
            if mod == n or (X, mod) in tried_complete or n // mod > COMPLETE_LIMIT:
                continue
            tried_complete.add((X, mod))
            a = complete_partial(oracle, X, r, mod, n)
            if a is not None:
                found[X] = a
                outcome[X] = "verified"
                log(f"  completed partial a = {r} mod {mod} -> verified {X} = g^{a}")
                if X in shot_strikes:
                    shot_strikes.remove(X)
            else:
                log(f"  {X}: no lift of a = {r} mod {mod} verifies")
                if X not in shot_strikes:
                    shot_strikes.append(X)

        if len(todo) >= 3 and all(outcome.get(X) in ("failed_full", "impossible") for X in todo):
            log(f"  joint shot {i}: every element failed -> blaming chi(g) (junk), no strikes")
            continue
        for X in shot_strikes:
            strikes[X] += 1
            log(f"  {X}: strike {strikes[X]}/{FORGERY_STRIKES}")
            if strikes[X] >= FORGERY_STRIKES:
                log(f"  {X} contradicts every exponent -> forged; stopping early")
                return found, X, n

    if not trimmed:                                       # never got a clean chi(g)
        n = reduce_to_order(oracle, "g", N)
    return found, None, n


# ------------------------------------------------------------------ entry

def solve(oracle):
    N = find_order(oracle, trim=False)
    if N is None:
        log("order finding failed; guessing 1")
        return {"order": 1, "secret": "none"}
    found, forged, order = find_exponents(oracle, N)
    if forged is None and all(X in found for X in ELEMENTS):
        secret = found["A"] * found["B"] * found["C"] % order
    else:
        if forged is None:
            missing = [X for X in ELEMENTS if X not in found]
            log(f"  never verified: {missing} -> treating as forgery")
        secret = "none"
    log(f"order = {order}, secret = {secret}, queries used = {oracle.queries_used}")
    return {"order": order, "secret": secret}