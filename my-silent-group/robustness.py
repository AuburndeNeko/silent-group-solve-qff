"""Robustness check: public instances + many mock groups, scored exactly like the judge.

    python my-silent-group/robustness.py my-silent-group/silent_group_solver_v2.py --seeds 3 --public 5
"""

import argparse
import importlib.util
import statistics
import time
from collections import defaultdict

from silentgroup import MockGroup, Oracle, load_public
from silentgroup.check import score
from silentgroup.oracle import BUDGET

SHAPES = [
    ("cyclic ~2^20", [1000003]),
    ("cyclic 2^59", [2**59]),
    ("cyclic ~2^58", [6145295, 53834531687]),
    ("non-cyclic small", [2**20, 6]),
    ("non-cyclic ~2^58", [6145295, 6145295 * 8761]),
    ("non-cyclic 3-power ~2^60", [3**18, 3**18 * 7]),
]
FORGERIES = [None, "C", "AB"]


def load_solver(path):
    spec = importlib.util.spec_from_file_location("submission", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.VERBOSE = False
    return module


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("solver")
    ap.add_argument("--seeds", type=int, default=3, help="mock seeds per shape/forgery")
    ap.add_argument("--public", type=int, default=5, help="runs per public instance")
    args = ap.parse_args()
    solver = load_solver(args.solver)

    results = defaultdict(list)
    failures = []

    def trial(label, instance, answer):
        oracle = Oracle(instance, budget=BUDGET)
        t0 = time.monotonic()
        try:
            res, err = solver.solve(oracle), None
        except Exception as e:
            res, err = None, repr(e)
        secs = time.monotonic() - t0
        pts = score(res, answer, BUDGET, oracle.queries_used)
        results[label].append((pts >= 1, pts, oracle.queries_used, secs))
        if pts < 1:
            failures.append((label, res, answer, err))
        print(f"{label:40s} {pts:.3f} pts  {oracle.queries_used:6d} q  {secs:5.1f}s", flush=True)

    for name in ("public_1", "public_2"):
        inst = load_public(name)
        for _ in range(args.public):
            trial(name, inst, inst.answer)

    for shape, moduli in SHAPES:
        for forged in FORGERIES:
            for seed in range(args.seeds):
                kwargs = {"seed": seed}
                if forged:
                    kwargs["forged"] = forged
                try:
                    mock = MockGroup(moduli, **kwargs)
                except ValueError:
                    continue            # forgery impossible for this group/seed
                label = f"{shape}, {'forged ' + forged if forged else 'genuine'}"
                trial(label, mock, mock.answer)

    print("\n" + "=" * 84)
    print(f"{'category':40s} {'correct':>9s} {'mean pts':>9s} {'mean q':>8s} {'max s':>7s}")
    all_rows = []
    for label, rows in results.items():
        all_rows += rows
        ok = f"{sum(r[0] for r in rows)}/{len(rows)}"
        print(f"{label:40s} {ok:>9s} {statistics.mean(r[1] for r in rows):9.3f} "
              f"{statistics.mean(r[2] for r in rows):8.0f} {max(r[3] for r in rows):7.1f}")
    n_ok = sum(r[0] for r in all_rows)
    print(f"\nOVERALL: {n_ok}/{len(all_rows)} correct, "
          f"mean {statistics.mean(r[1] for r in all_rows):.3f} pts, "
          f"mean {statistics.mean(r[2] for r in all_rows):.0f} queries, "
          f"max {max(r[3] for r in all_rows):.1f}s")
    for label, res, answer, err in failures:
        print(f"\nFAILED {label}: got {res}, expected {answer}" + (f", error {err}" if err else ""))


if __name__ == "__main__":
    main()