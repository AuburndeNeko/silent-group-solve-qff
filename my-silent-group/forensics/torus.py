"""Character torus in 3-D: genuine broadcasts lie on a torus knot, forged ones don't.

Each shot samples one character chi. (chi(g), chi(X)) are two angles, i.e. a point on a
torus: chi(g) goes around the big ring, chi(X) around the tube. For genuine X = g^a,
chi(X) = a*chi(g) mod 1, so the points lie on a (1, a) torus knot. Saves torus3d.png
(and torus3d.gif with --gif).

    python my-silent-group/torus3d.py my-silent-group/submission.py [--gif]
"""

import importlib.util
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from silentgroup import MockGroup, Oracle

spec = importlib.util.spec_from_file_location("solver", sys.argv[1])
S = importlib.util.module_from_spec(spec)
spec.loader.exec_module(S)
S.VERBOSE = False

M, SHOTS, BITS = 31, 300, 14
R, r = 2.0, 0.8                      # ring radius, tube radius

mock = None
for seed in range(50):
    try:
        mock = MockGroup([M, M], seed=seed, forged="C")
        break
    except ValueError:
        continue
if mock is None:
    sys.exit("could not build a forged mock; try other moduli")

oracle = Oracle(mock, noise=0)
names = ["g", "A", "C"]
qc = S.ipe_circuit(names, BITS)
pts = {"A": [], "C": []}
for _ in range(SHOTS):
    ph = S.read_phases(oracle.run(qc, shots=1), names, BITS)
    for X in pts:
        pts[X].append((float(ph["g"]), float(ph[X])))


def embed(u, v):
    """(u, v) in [0,1)^2 -> point on a torus in 3-D."""
    th, phi = 2 * np.pi * np.asarray(u), 2 * np.pi * np.asarray(v)
    return ((R + r * np.cos(phi)) * np.cos(th),
            (R + r * np.cos(phi)) * np.sin(th),
            r * np.sin(phi))


U, V = np.meshgrid(np.linspace(0, 1, 72), np.linspace(0, 1, 36))
fig = plt.figure(figsize=(12, 6))
axes = []
for i, X in enumerate(["A", "C"]):
    ax = fig.add_subplot(1, 2, i + 1, projection="3d")
    axes.append(ax)
    ax.plot_surface(*embed(U, V), color="lightgray", alpha=0.12, linewidth=0)
    ax.plot_wireframe(*embed(U, V), color="gray", alpha=0.07, rstride=3, cstride=3)
    xs, ys = zip(*pts[X])
    grid = {(round(x * M) % M, round(y * M) % M) for x, y in pts[X]}
    if X == "A":
        xg, xa = next((u, v) for u, v in grid if u != 0)
        a = xa * pow(xg, -1, M) % M
        t = np.linspace(0, 1, 4000)
        ax.plot(*embed(t, (a * t) % 1), lw=1.0, alpha=0.7, color="tab:orange")
        title = f"A = g^{a} (genuine): {len(grid)} points on a (1, {a}) torus knot"
    else:
        title = f"C forged: {len(grid)} points scattered over the surface"
    ax.scatter(*embed(xs, ys), s=16, depthshade=False, color="tab:blue")
    ax.set_title(title, fontsize=10)
    ax.set_box_aspect((1, 1, 0.45))
    ax.set_axis_off()
    ax.view_init(elev=38, azim=30)

fig.suptitle(f"Character torus: mock group Z_{M} x Z_{M}, {SHOTS} shots, no noise")
fig.tight_layout()
fig.savefig("torus3d.png", dpi=150)
print("saved torus3d.png")

if "--gif" in sys.argv:
    from matplotlib.animation import FuncAnimation, PillowWriter

    def turn(frame):
        for ax in axes:
            ax.view_init(elev=38, azim=30 + 4 * frame)
        return axes

    FuncAnimation(fig, turn, frames=90, interval=60).save("torus3d.gif", writer=PillowWriter(fps=15))
    print("saved torus3d.gif")