"""Open track: forgery forensics + partial key recovery.

Per broadcast X in A, B, C:
  genuine -> its exponent (X = g^a), verified on the device
  forged  -> ord(X); smallest q with X^q in <g>; verified X^q = g^w; <g,X> = Z_d1 x Z_d2;
             split X = g^c * y with y^q = e (when it exists)
Across forged broadcasts: the tag group Q = <g,A,B,C>/<g> (order, exponent, verified relations
  between tags, e.g. A*B in <g> means the same tag was reused). Combined elements need no new
  circuits: chi(A^i B^j) = i chi(A) + j chi(B) is computed from the same joint shots.
About the hidden group G: |G| >= |<g,A,B,C>| = n|Q|, and a certificate that G is non-cyclic.
Partial key recovery: pair keys of genuine parties, verified on the device.
No query budget: this is analysis, not the scored solver.

    python my-silent-group/forensics.py my-silent-group/submission.py [--repeat K]
"""

import importlib.util
import itertools
import math
import sys
from fractions import Fraction

from silentgroup import MockGroup, Oracle, load_public

GREEN, RED, YELLOW, CYAN, BOLD, DIM, RESET = (
    "\033[32m", "\033[31m", "\033[33m", "\033[36m", "\033[1m", "\033[2m", "\033[0m")

spec = importlib.util.spec_from_file_location("solver", sys.argv[1])
S = importlib.util.module_from_spec(spec)
spec.loader.exec_module(S)
S.VERBOSE = False

JOINT_SHOTS = 10
LIFT_LIMIT = 64
TAG_BOX_CAP = 2000       # max size of the search box for tag relations
TAG_VERIFY_CAP = 40      # max device verifications in the tag search


def divisors(n):
    ds = [1]
    for p, e in S.factorize(n).items():
        ds = [d * p ** k for d in ds for k in range(e + 1)]
    return sorted(ds)


def joint_readings(oracle, names, D, shots):
    """Per shot: {name: x} with chi(name) = x / D; clearly-junk shots dropped."""
    m = D.bit_length() + S.MARGIN_BITS
    qc = S.ipe_circuit(names, m)
    out = []
    for _ in range(shots):
        ph = S.read_phases(oracle.run(qc, shots=1), names, m)
        ys = {nm: ph[nm] * D for nm in names}
        if all(abs(y - round(y)) <= Fraction(1, 4) for y in ys.values()):
            out.append({nm: round(y) % D for nm, y in ys.items()})
    return out


def merged_congruence(expo, readings, D):
    """Solve w*x_g = sum e_X x_X (mod D) across shots (one junk shot tolerated).
    Returns (r, mod) meaning w = r (mod mod), or None."""
    congs = [S.congruence(rd["g"], sum(e * rd[X] for X, e in expo.items()) % D, D)
             for rd in readings]
    for sub in [congs] + [congs[:i] + congs[i + 1:] for i in range(len(congs))]:
        if not sub or any(c is None for c in sub):
            continue
        acc = (0, 1)
        for c in sub:
            acc = S.crt_merge(acc, c)
            if acc is None:
                break
        if acc is not None:
            return acc
    return None


def find_w(oracle, expo, readings, n, D):
    """w with prod X^e = g^w, verified on the device; None if not found."""
    acc = merged_congruence(expo, readings, D)
    if acc is None:
        return None
    r, mod = acc
    if n // mod > LIFT_LIMIT:
        return None
    for j in range(n // mod):
        w = (r + j * mod) % n
        terms = [(X, e) for X, e in expo.items() if e] + ([("g", -w)] if w else [])
        if (S.is_identity(oracle, terms, max_shots=8, batch=4, reject_at=1)
                and S.is_identity(oracle, terms)):
            return w
    return None


def tag_group(oracle, F, q, readings, n, D):
    """Q = <g, forged>/<g>: (|Q|, exponent of Q, [(expo, w)] verified relations) or None."""
    size = math.prod(q[X] for X in F)
    if size > TAG_BOX_CAP:
        return None
    kernel, relations, checks = {tuple(0 for _ in F)}, [], 0
    for e in itertools.product(*[range(q[X]) for X in F]):
        if not any(e):
            continue
        expo = dict(zip(F, e))
        if merged_congruence(expo, readings, D) is None:
            continue                                  # classically impossible: skip device
        if checks >= TAG_VERIFY_CAP:
            break
        checks += 1
        w = find_w(oracle, expo, readings, n, D)
        if w is not None:
            kernel.add(e)
            relations.append((expo, w))
    order_Q = size // len(kernel)

    def elt_order(e):
        k = 1
        while tuple(k * ei % q[X] for ei, X in zip(e, F)) not in kernel:
            k += 1
        return k

    exponent = max(elt_order(e) for e in itertools.product(*[range(q[X]) for X in F]))
    closed = all(tuple((a + b) % q[X] for a, b, X in zip(e1, e2, F)) in kernel
                 for e1 in kernel for e2 in kernel)
    return order_Q, exponent, relations, closed


def fmt_expo(expo):
    return " * ".join(f"{X}^{e}" if e != 1 else X for X, e in expo.items() if e)


def investigate(label, inst, moduli=None, quiet=False):
    say = (lambda *a: None) if quiet else print
    oracle = Oracle(inst)
    n = S.find_order(oracle, "g", trim=True)
    say(f"\n=== {label} ===")
    say(f"ord(g) = n = {n}")
    if moduli:
        a, b = moduli
        say(f"(ground truth: G = Z_{math.gcd(a, b)} x Z_{a * b // math.gcd(a, b)})")
    summary = {"n": n}
    exps, forged, q, noncyclic = {}, [], {}, []

    for X in S.ELEMENTS:
        found, _, _ = S.find_exponents(oracle, n, elements=(X,))
        if X in found:
            exps[X] = found[X]
            summary[X] = ("genuine", found[X])
            say(f"  {BOLD}{X}{RESET} = g^{found[X]}   {GREEN}(genuine, verified on device){RESET}")
            continue
        nX = S.find_order(oracle, X, trim=True)
        D = S.lcm(n, nX)
        readings = joint_readings(oracle, ["g", X], D, JOINT_SHOTS)
        qX = w = None
        for cand in divisors(nX):
            w = find_w(oracle, {X: cand}, readings, n, D)
            if w is not None:
                qX = cand
                break
        if qX is None:
            summary[X] = ("forged", nX, None)
            say(f"  {BOLD}{X}{RESET}: {RED}{BOLD}FORGED{RESET}; ord({X}) = {nX}; "
                f"{DIM}relation to <g> not determined{RESET}")
            continue
        if qX == 1:
            exps[X] = w
            summary[X] = ("genuine", w)
            say(f"  {BOLD}{X}{RESET} = g^{w}   {GREEN}(genuine after all, verified on device){RESET}")
            continue
        forged.append(X)
        q[X] = qX
        d1 = math.gcd(math.gcd(n, w), qX)
        d2 = n * qX // d1
        if d1 > 1:
            noncyclic.append((X, d1, d2))
        gq = math.gcd(qX, n)
        split = None
        if w % gq == 0:
            split = ((w // gq) * pow(qX // gq, -1, n // gq) % (n // gq), n // gq)
        summary[X] = ("forged", nX, qX, w, d1, d2, split)
        say(f"  {BOLD}{X}{RESET}: {RED}{BOLD}FORGED{RESET} (not a power of g)")
        say(f"     {DIM}ord({X}) ={RESET} {nX}")
        say(f"     smallest q with {X}^q in <g>: q = {YELLOW}{BOLD}{qX}{RESET};  "
            f"{GREEN}verified{RESET} {X}^{qX} = g^{YELLOW}{w}{RESET}")
        say(f"     [<g,{X}> : <g>] = {YELLOW}{qX}{RESET};   "
            f"<g,{X}> = {CYAN}Z_{d1} x Z_{d2}{RESET}  {DIM}(order {n * qX}){RESET}")
        if moduli:
            a, b = moduli
            G1, G2 = math.gcd(a, b), a * b // math.gcd(a, b)
            fits = G1 % d1 == 0 and G2 % d2 == 0
            say(f"     fits inside G = {CYAN}Z_{G1} x Z_{G2}{RESET}: "
                f"{GREEN + 'yes' if fits else RED + 'NO'}{RESET}")
        if split:
            say(f"     split: {X} = g^{YELLOW}{split[0]}{RESET} * y  with  y^{qX} = e   "
                f"{DIM}(c defined mod {split[1]}){RESET}")
        else:
            say(f"     {DIM}no split X = g^c * y with y^{qX} = e   "
                f"(gcd(q, n) = {gq} does not divide w){RESET}")

    order_Q = None
    if len(forged) >= 2:
        nF = [S.find_order(oracle, X, trim=True) for X in forged]
        D = n
        for v in nF:
            D = S.lcm(D, v)
        readings = joint_readings(oracle, ["g"] + forged, D, JOINT_SHOTS)
        tg = tag_group(oracle, forged, q, readings, n, D)
        say(f"  tag group Q = <g,{','.join(forged)}>/<g>:")
        if tg is None:
            say(f"     search box too large (> {TAG_BOX_CAP}); skipped")
        else:
            order_Q, expQ, rels, closed = tg
            summary["Q"] = (order_Q, expQ, tuple(sorted((tuple(r[0].items()), r[1]) for r in rels)))
            shape = "cyclic" if expQ == order_Q else f"non-cyclic (exponent {expQ})"
            say(f"     |Q| = {order_Q}, {shape}; product of tag orders = {math.prod(q.values())}")
            for expo, w in rels:
                say(f"     verified relation: {fmt_expo(expo)} = g^{w}   (in <g>: tags cancel)")
            if not rels:
                say("     no relations: the tags are independent (different fakes)")
            if not closed:
                say("     warning: found relations are not closed under addition (a check failed?)")
    elif len(forged) == 1:
        order_Q = q[forged[0]]

    say("  about the hidden group G:")
    if order_Q is not None:
        say(f"     |G| >= |<g,A,B,C>| = n * |Q| = {n * order_Q}  (~2^{math.log2(n * order_Q):.1f})")
    else:
        say(f"     |G| >= n = {n}")
    if noncyclic:
        X, d1, d2 = noncyclic[0]
        say(f"     G is NOT cyclic: <g,{X}> = Z_{d1} x Z_{d2} with gcd {d1} > 1, and every "
            f"subgroup of a cyclic group is cyclic")
    else:
        say("     no certificate of non-cyclicity from these broadcasts")

    say("  partial key recovery:")
    for X, Y in itertools.combinations(sorted(exps), 2):
        k = exps[X] * exps[Y] % n
        terms = [(X, exps[Y])] + ([("g", -k)] if k else [])
        ok = S.is_identity(oracle, terms)
        summary[f"key_{X}{Y}"] = k
        say(f"     pair key g^({X.lower()}{Y.lower()}) = g^{k}   "
            f"device check {X}^{Y.lower()} = g^({X.lower()}{Y.lower()}): {'OK' if ok else 'FAILED'}")
    if len(exps) == 3:
        say(f"     full secret abc mod n = {exps['A'] * exps['B'] * exps['C'] % n}")
    if not exps:
        say("     no genuine parties")
    say(f"  (queries used for this analysis: {oracle.queries_used})")
    return summary


CASES = [
    ("public_2", lambda: load_public("public_2"), None),
    ("mock Z_2^20 x Z_6, forged AB", lambda: MockGroup([2**20, 6], seed=1, forged="AB"), (2**20, 6)),
    ("mock non-cyclic ~2^58, forged C",
     lambda: MockGroup([6145295, 6145295 * 8761], seed=2, forged="C"), (6145295, 6145295 * 8761)),
    ("public_1 (all genuine: control)", lambda: load_public("public_1"), None),
]

if __name__ == "__main__":
    if "--repeat" in sys.argv:
        K = int(sys.argv[sys.argv.index("--repeat") + 1])
        print(f"reproducibility: public_2 analysed {K} times")
        first = investigate("public_2", load_public("public_2"))
        same = 1
        for i in range(1, K):
            s = investigate("public_2", load_public("public_2"), quiet=True)
            same += s == first
            print(f"  run {i + 1}: {'identical' if s == first else 'DIFFERENT: ' + str(s)}")
        print(f"{same}/{K} runs gave identical results")
    else:
        for label, make, moduli in CASES:
            investigate(label, make(), moduli)