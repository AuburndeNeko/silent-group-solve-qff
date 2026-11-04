"""Silent Group solver v3 (experimental).

v2.1, except the first shot is a "fat shot": it reads g, A, B, C at FAT_BITS each and a
small lattice (LLL) recovers their common denominator, tolerating one junk reading. One
clean fat shot gives an order candidate AND exponent congruences, replacing the 121-bit
order readout plus the first joint shot. Anything unsettled falls back to v2.1 steps.
"""

import math
import random
from fractions import Fraction

from qiskit import ClassicalRegister, QuantumCircuit, QuantumRegister

from silentgroup import oracle_gate

MAX_BITS = 60        # |G| <= 2^60, so every element's order is <= 2^60
MARGIN_BITS = 4      # extra bits when rounding to the 1/n grid
FORGERY_STRIKES = 3  # contradicting shots before calling an element forged
FAT_BITS = 82        # bits per element in the fat shot
LATTICE_W = 32        # weight on the denominator coordinate of the lattice
ERR_TOL = 8         # max readout error (units of 2^-FAT_BITS) for a reading to count as clean
VERBOSE = True       # set False before submitting
ELEMENTS = ("A", "B", "C")


def log(*args):
    if VERBOSE:
        print(*args)


def can_afford(oracle, qc, shots=1):
    left = oracle.budget_left()
    return left is None or oracle.cost(qc, shots) <= left


# ---------------------------------------------------------------- circuits

def add_ipe(qc, qubit, cbits, element):
    """Iterative phase estimation of chi(element) to len(cbits) bits; b_k -> cbits[k-1]."""
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
    q = QuantumRegister(1, "q")
    c = ClassicalRegister(m * len(elements), "c")
    qc = QuantumCircuit(q, c)
    for j, el in enumerate(elements):
        add_ipe(qc, q[0], [c[j * m + i] for i in range(m)], el)
    return qc


def read_values(counts, elements, m):
    """1-shot counts -> {element: integer v}; the phase is exactly v / 2^m."""
    (key,) = counts
    bits = key[::-1]
    return {el: int(bits[j * m:(j + 1) * m], 2) for j, el in enumerate(elements)}


def read_phases(counts, elements, m):
    return {el: Fraction(v, 2 ** m) for el, v in read_values(counts, elements, m).items()}


def identity_test_circuit(terms):
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
    """From xX = a*xg (mod n): (r, mod) meaning a = r (mod mod), or None if impossible."""
    d = math.gcd(xg, n)
    if xX % d:
        return None
    mod = n // d
    if mod == 1:
        return (0, 1)
    return ((xX // d) * pow(xg // d, -1, mod) % mod, mod)


def crt_merge(c1, c2):
    (r1, m1), (r2, m2) = c1, c2
    g = math.gcd(m1, m2)
    if (r2 - r1) % g:
        return None
    l = m1 // g * m2
    m2g = m2 // g
    t = 0 if m2g == 1 else ((r2 - r1) // g) * pow(m1 // g, -1, m2g) % m2g
    return ((r1 + m1 * t) % l, l)


# ---------------------------------------------------------------- lattice

def lll_reduce(basis, delta=Fraction(99, 100)):
    """Textbook LLL with exact rational Gram-Schmidt (fine for our 4x4 / 5x5 bases)."""
    b = [list(row) for row in basis]
    n = len(b)

    def gso():
        bs, nrm = [], []
        mu = [[Fraction(0)] * n for _ in range(n)]
        for i in range(n):
            v = [Fraction(x) for x in b[i]]
            for j in range(i):
                mu[i][j] = sum(Fraction(x) * y for x, y in zip(b[i], bs[j])) / nrm[j]
                v = [vi - mu[i][j] * yj for vi, yj in zip(v, bs[j])]
            bs.append(v)
            nrm.append(sum(x * x for x in v))
        return mu, nrm

    mu, nrm = gso()
    k = 1
    while k < n:
        for j in range(k - 1, -1, -1):          # size-reduce b[k]
            q = round(mu[k][j])
            if q:
                b[k] = [x - q * y for x, y in zip(b[k], b[j])]
                for i in range(j):
                    mu[k][i] -= q * mu[j][i]
                mu[k][j] -= q
        if nrm[k] >= (delta - mu[k][k - 1] ** 2) * nrm[k - 1]:
            k += 1
        else:                                   # Lovasz condition fails: swap
            b[k], b[k - 1] = b[k - 1], b[k]
            mu, nrm = gso()
            k = max(k - 1, 1)
    return b


def lattice_candidates(vals, M):
    """Candidate common denominators q <= 2^MAX_BITS for phases vals[i] / 2^M.

    Lattice rows: (W, v1, ..., vk) and 2^M * e_i. The vector q*row0 - sum p_i*row_i is
    (qW, q*v1 - p1*2^M, ...): short exactly when q makes every q*phase nearly whole.
    """
    k = len(vals)
    basis = [[LATTICE_W] + list(vals)]
    for i in range(k):
        row = [0] * (k + 1)
        row[i + 1] = 2 ** M
        basis.append(row)
    red = lll_reduce(basis)
    vecs = list(red)
    for i in range(len(red)):
        for j in range(i + 1, len(red)):
            vecs.append([x + y for x, y in zip(red[i], red[j])])
            vecs.append([x - y for x, y in zip(red[i], red[j])])
    out = set()
    for v in vecs:
        q = abs(v[0]) // LATTICE_W
        if 0 < q <= 2 ** MAX_BITS:
            out.add(q)
    return out


def fits(q, v, M):
    """Is q * (v / 2^M) within q*ERR_TOL / 2^M of a whole number?"""
    res = (q * v) % (2 ** M)
    return min(res, 2 ** M - res) <= q * ERR_TOL


def fat_shot_denominator(vals, M):
    """Common denominator of the readings, tolerating one bad reading.

    Tries the full set and every 3-subset; keeps the q that fits the MOST readings
    (ties -> smallest q). Returns (q, clean_names) or (None, []).
    """
    names = list(vals)
    subsets = [names] + [[x for x in names if x != drop] for drop in names]
    best = None
    for sub in subsets:
        for q in lattice_candidates([vals[x] for x in sub], M):
            if not all(fits(q, vals[x], M) for x in sub):
                continue
            clean = [x for x in names if fits(q, vals[x], M)]
            key = (len(clean), -q)
            if best is None or key > best[0]:
                best = (key, q, clean)
        if sub is names and best is not None and len(best[2]) == len(names):
            break                       # all readings clean: no need to try subsets
    return (best[1], best[2]) if best else (None, [])


# ---------------------------------------------------- quantum subroutines

def is_identity(oracle, terms, max_shots=16, batch=4):
    qc = identity_test_circuit(terms)
    ones = shots = 0
    while shots < max_shots:
        counts = oracle.run(qc, shots=batch)
        ones += counts.get("1", 0)
        shots += batch
        if ones >= 2:
            return False
    return True


def verify(oracle, X, a):
    """Is X == g^a?"""
    return is_identity(oracle, [(X, 1)] + ([("g", -a)] if a else []))


def reduce_to_order(oracle, element, N):
    for p in sorted(factorize(N)):
        while N % p == 0 and is_identity(oracle, [(element, N // p)]):
            N //= p
    return N


class OrderSearch:
    """Keeps every lcm of subsets of denominator readings; tests them with g^N = e."""

    def __init__(self, oracle, element="g"):
        self.oracle, self.element = oracle, element
        self.candidates, self.failed = set(), []

    def add(self, r):
        new = {r} | {lcm(C, r) for C in self.candidates}
        new = {N for N in new if N <= 2 ** MAX_BITS} - self.candidates
        self.candidates |= new
        for N in sorted(new, reverse=True):
            if any(F % N == 0 for F in self.failed):
                continue
            if is_identity(self.oracle, [(self.element, N)]):
                order = reduce_to_order(self.oracle, self.element, N)
                log(f"  g^{N} = e passed; after trimming, order = {order}")
                return order
            self.failed.append(N)
        return None


def find_order(oracle, element="g", max_shots=8, search=None):
    search = search or OrderSearch(oracle, element)
    m = 2 * MAX_BITS + 1
    qc = ipe_circuit([element], m)
    for i in range(max_shots):
        if not can_afford(oracle, qc):
            log("  out of budget for order finding")
            break
        phi = read_phases(oracle.run(qc, shots=1), [element], m)[element]
        r = phi.limit_denominator(2 ** MAX_BITS).denominator
        log(f"  121-bit shot {i}: r = {r}")
        order = search.add(r)
        if order:
            return order
    return None


def find_exponents(oracle, n, elements=ELEMENTS, max_shots=6, init_cands=None, init_strikes=None):
    """v2.1 exponent search, optionally seeded with congruences/strikes from the fat shot."""
    m = n.bit_length() + MARGIN_BITS
    cands = {X: set((init_cands or {}).get(X, ())) for X in elements}
    failed = {X: set() for X in elements}
    strikes = {X: (init_strikes or {}).get(X, 0) for X in elements}
    found = {}

    for X in elements:                          # verify what the fat shot already pinned down
        for r in sorted(r for r, mod in cands[X] if mod == n):
            if verify(oracle, X, r):
                found[X] = r
                log(f"  verified from fat shot: {X} = g^{r}")
                break
            failed[X].add(r)
        if X not in found and failed[X]:
            strikes[X] += 1
    for X in elements:
        if X not in found and strikes[X] >= FORGERY_STRIKES:
            return found, X

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
                    if verify(oracle, X, r):
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
    names = ["g"] + list(ELEMENTS)
    search = OrderSearch(oracle)
    order, fat = None, None
    strikes0 = {X: 0 for X in ELEMENTS}

    qc = ipe_circuit(names, FAT_BITS)
    if can_afford(oracle, qc):
        vals = read_values(oracle.run(qc, shots=1), names, FAT_BITS)
        q, clean = fat_shot_denominator(vals, FAT_BITS)
        log(f"  fat shot: common denominator q = {q}, clean = {clean}")
        if q is not None:
            fat = (q, {x: round(Fraction(vals[x] * q, 2 ** FAT_BITS)) % q for x in clean})
            if "g" in clean:
                for X in ELEMENTS:
                    if X not in clean:
                        strikes0[X] += 1
            order = search.add(q)

    if order is None:
        log("  order not settled by the fat shot -> 121-bit readouts")
        order = find_order(oracle, search=search)
    if order is None:
        log("order finding failed; guessing 1")
        return {"order": 1, "secret": "none"}

    init = {X: set() for X in ELEMENTS}
    if fat is not None:
        q, xs = fat
        if order % q == 0 and "g" in xs:
            for X in ELEMENTS:
                if X in xs:
                    cong = congruence(xs["g"], xs[X], q)
                    if cong is None:
                        strikes0[X] += 1
                    else:
                        init[X].add(cong)
                        log(f"  fat shot: {X}: a = {cong[0]} mod {cong[1]}")

    found, forged = find_exponents(oracle, order, init_cands=init, init_strikes=strikes0)
    if forged is None and all(X in found for X in ELEMENTS):
        secret = found["A"] * found["B"] * found["C"] % order
    else:
        if forged is None:
            missing = [X for X in ELEMENTS if X not in found]
            log(f"  never verified: {missing} -> treating as forgery")
        secret = "none"
    log(f"order = {order}, secret = {secret}, queries used = {oracle.queries_used}")
    return {"order": order, "secret": secret}