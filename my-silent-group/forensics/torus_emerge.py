"""Watch the forger get caught: points arrive shot by shot on the character torus.

Both panels draw the knot predicted from the first shot (exponent guessed from it).
Green = lands on the predicted knot, red = misses. Genuine: all green. Forged: red piles up.
Saves torus_emerge.gif and torus_emerge_final.png.

    python my-silent-group/torus_emerge.py my-silent-group/submission.py
"""

import importlib.util
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter

from silentgroup import MockGroup, Oracle

spec = importlib.util.spec_from_file_location("solver", sys.argv[1])
S = importlib.util.module_from_spec(spec)
spec.loader.exec_module(S)
S.VERBOSE = False

M, SHOTS, BITS = 31, 120, 14
PER_FRAME, HOLD = 2, 20              # points added per frame; frames held at the end
R, r = 2.0, 0.8

mock = None
for seed in range(50):
    try:
        mock = MockGroup([M, M], seed=seed, forged="C")
        break
    except ValueError:
        continue
if mock is None:
    sys.exit("could not build a forged mock")

oracle = Oracle(mock, noise=0)
names = ["g", "A", "C"]
qc = S.ipe_circuit(names, BITS)
shots = []
for _ in range(SHOTS):
    ph = S.read_phases(oracle.run(qc, shots=1), names, BITS)
    shots.append({k: round(v * M) % M for k, v in ph.items()})


def guess_exponent(X):
    for s in shots:
        if s["g"]:
            a = s[X] * pow(s["g"], -1, M) % M
            return a - M if a > M // 2 else a          # fewest windings
    return 0


def embed(u, v):
    th, phi = 2 * np.pi * np.asarray(u, dtype=float), 2 * np.pi * np.asarray(v, dtype=float)
    return ((R + r * np.cos(phi)) * np.cos(th),
            (R + r * np.cos(phi)) * np.sin(th),
            r * np.sin(phi))


U, V = np.meshgrid(np.linspace(0, 1, 72), np.linspace(0, 1, 36))
t = np.linspace(0, 1, 4000)
fig = plt.figure(figsize=(12, 6))
panels = []
for i, (X, label) in enumerate([("A", "A (genuine)"), ("C", "C (forged)")]):
    ax = fig.add_subplot(1, 2, i + 1, projection="3d")
    a = guess_exponent(X)
    ax.plot_surface(*embed(U, V), color="lightgray", alpha=0.12, linewidth=0)
    ax.plot(*embed(t, (a * t) % 1), lw=1.0, alpha=0.6, color="tab:orange")
    sc = ax.scatter([], [], [], s=18, depthshade=False)
    txt = ax.text2D(0.02, 0.95, "", transform=ax.transAxes, fontsize=10)
    ax.set_title(f"{label}: knot predicted from shot 1 (exponent {a})", fontsize=10)
    ax.set_box_aspect((1, 1, 0.45))
    ax.set_axis_off()
    panels.append((X, a, ax, sc, txt))
fig.suptitle("Catching a forger: green = on the predicted knot, red = off it")


def update(frame):
    k = min(SHOTS, (frame + 1) * PER_FRAME)
    sub = shots[:k]
    for X, a, ax, sc, txt in panels:
        xs, ys, zs = embed([s["g"] / M for s in sub], [s[X] / M for s in sub])
        sc._offsets3d = (xs, ys, zs)
        on = [(s[X] - a * s["g"]) % M == 0 for s in sub]
        sc.set_color(["tab:green" if o else "tab:red" for o in on])
        txt.set_text(f"shots: {k}    on the knot: {sum(on)}/{k}")
        ax.view_init(elev=38, azim=30 + 2 * frame)
    return [p[3] for p in panels]


frames = SHOTS // PER_FRAME + HOLD
FuncAnimation(fig, update, frames=frames, interval=80).save(
    "torus_emerge.gif", writer=PillowWriter(fps=12))
update(frames - 1)
fig.savefig("torus_emerge_final.png", dpi=150)
print("saved torus_emerge.gif and torus_emerge_final.png")