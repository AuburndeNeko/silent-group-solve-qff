"""Classical self-test of the fat-shot lattice: no device, no queries.

    python my-silent-group/lattice_selftest.py my-silent-group/silent_group_solver_v3.py
"""

import importlib.util
import math
import random
import sys
import time

spec = importlib.util.spec_from_file_location("solver", sys.argv[1])
solver = importlib.util.module_from_spec(spec)
spec.loader.exec_module(solver)
solver.VERBOSE = False

ERRORS = [0] * 6 + [1, -1] * 3 + [2, -2, 3, -3]   # readout error, in last-bit units
TRIALS = 50


def make_shot(log_r, M, junk):
    r = random.randrange(2 ** (log_r - 1), 2 ** log_r)
    while True:
        xg = random.randrange(1, r)
        if math.gcd(xg, r) == 1:
            break
    nums = {"g": xg}
    for X in "ABC":
        nums[X] = random.randrange(r) * xg % r            # chi(X) = a * chi(g)
    vals = {}
    for name, x in nums.items():
        v = (2 * x * 2 ** M + r) // (2 * r)                # nearest integer to x * 2^M / r
        vals[name] = (v + random.choice(ERRORS)) % 2 ** M
    if junk:
        vals[junk] = random.randrange(2 ** M)
    return r, vals


print(f"{'log2 r':>6} {'M':>4} {'junk':>5} {'success':>9} {'ms/shot':>8}")
for log_r in (40, 56, 58, 60):
    for M in (82, 84, 86):
        for junk in (None, "B", "g"):
            ok = 0
            t0 = time.monotonic()
            for _ in range(TRIALS):
                r, vals = make_shot(log_r, M, junk)
                q, clean = solver.fat_shot_denominator(vals, M)
                ok += (q == r) and (junk is None or junk not in clean)
            ms = (time.monotonic() - t0) / TRIALS * 1000
            print(f"{log_r:6d} {M:4d} {str(junk):>5} {ok:5d}/{TRIALS:<3d} {ms:8.1f}")