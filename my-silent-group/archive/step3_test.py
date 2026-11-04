import time

from silentgroup import MockGroup, Oracle

import silent_group_solver as solver

cases = [
    ("cyclic prime",           MockGroup([1000003], g=[5], seed=0)),
    ("non-cyclic small",       MockGroup([2**20, 6], seed=1)),
    ("public-sized",           MockGroup([6145295, 53834531687], seed=2)),
    ("non-cyclic big, C forged", MockGroup([6145295, 6145295 * 8761], seed=2, forged="C")),
    ("non-cyclic, AB forged", MockGroup([2 ** 20, 6], seed=1, forged="AB")),
    ("non-cyclic big, genuine", MockGroup([6145295, 6145295 * 8761], seed=2)),
]

for label, mock in cases:
    print(f"\n=== {label} ===")
    oracle = Oracle(mock)
    t0 = time.monotonic()
    res = solver.solve(oracle)
    ok = res["order"] == mock.answer["order"] and res["secret"] == mock.answer["secret"]
    print(f"{'OK' if ok else 'WRONG'}: got {res}, true {mock.answer}, "
          f"queries {oracle.queries_used}, {time.monotonic() - t0:.1f}s")