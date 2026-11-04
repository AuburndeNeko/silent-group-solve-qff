"""Forger fingerprints: barcodes of the prime factors of each forgery's tag order q.

Band width = exponent * log2(p) (its size in bits); each prime keeps one colour across all
forgeries, so equal fingerprints stand out. A thinner bar underneath shows ord(X).
Values come from forensics.py runs. Saves fingerprint.png.

    python my-silent-group/fingerprint.py my-silent-group/submission.py
"""

import importlib.util
import math
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

spec = importlib.util.spec_from_file_location("solver", sys.argv[1])
S = importlib.util.module_from_spec(spec)
spec.loader.exec_module(S)

# (label, ord(X), tag order q) from the forensics output
CASES = [
    ("public_2: B", 1771394625, 585585),
    ("mock Z_2^20 x Z_6: A", 32768, 2),
    ("mock Z_2^20 x Z_6: B", 262144, 2),
    ("mock ~2^58: C", 53838929495, 6145295),
]

facs = [(lab, S.factorize(nX), S.factorize(q)) for lab, nX, q in CASES]
primes = sorted({p for _, fx, fq in facs for p in list(fx) + list(fq)})
cmap = plt.get_cmap("tab20")
colour = {p: cmap(i % 20) for i, p in enumerate(primes)}

fig, ax = plt.subplots(figsize=(11, 1.4 * len(CASES) + 1))
for row, (lab, fx, fq) in enumerate(facs):
    y = len(CASES) - 1 - row
    for fac, dy, h, big in ((fq, 0.12, 0.42, True), (fx, -0.28, 0.16, False)):
        x = 0.0
        for p in sorted(fac):
            w = fac[p] * math.log2(p)
            ax.barh(y + dy, w, left=x, height=h, color=colour[p], edgecolor="white")
            if big and w > 1.6:
                txt = f"{p}" + (f"^{fac[p]}" if fac[p] > 1 else "")
                ax.text(x + w / 2, y + dy, txt, ha="center", va="center", fontsize=9)
            x += w
    q_val = math.prod(p ** e for p, e in fq.items())
    ax.text(-0.5, y + 0.12, f"{lab}\nq = {q_val}", ha="right", va="center", fontsize=9)

ax.set_yticks([])
ax.set_xlabel("bits (exponent x log2 p)    thick bar: tag order q    thin bar: ord(X)")
ax.set_xlim(left=0)
ax.set_ylim(-0.6, len(CASES) - 0.3)
for side in ("top", "right", "left"):
    ax.spines[side].set_visible(False)
ax.set_title("Forger fingerprints: prime factorization of each forgery's tag order")
fig.tight_layout()
fig.savefig("fingerprint.png", dpi=150, bbox_inches="tight")
print("saved fingerprint.png")
for lab, fx, fq in facs:
    print(f"  {lab}: q = " + " * ".join(f"{p}^{e}" if e > 1 else str(p) for p, e in sorted(fq.items())))

PALETTE = [196, 208, 226, 46, 51, 27, 129, 201, 214, 34, 93, 160]   # 256-colour codes
ansi = {p: PALETTE[i % len(PALETTE)] for i, p in enumerate(primes)}


def bar(fac):
    out = ""
    for p in sorted(fac):
        w = max(1, round(fac[p] * math.log2(p)))
        label = f"{p}" + (f"^{fac[p]}" if fac[p] > 1 else "")
        cell = label.center(w) if w >= len(label) else " " * w
        out += f"\033[48;5;{ansi[p]}m\033[38;5;16m{cell}\033[0m"
    return out


print()
for lab, fx, fq in facs:
    print(f"\033[1m{lab:24s}\033[0m q      {bar(fq)}")
    print(f"{'':24s} ord(X) {bar(fx)}")
    print()