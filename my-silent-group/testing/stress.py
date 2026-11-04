"""Run one instance many times; print the full log only for failed runs.

    python my-silent-group/stress.py my-silent-group/silent_group_solver_v4.1.py pow2 100

Instances: public_1, public_2, or a mock: pow2, pow3.
"""

import contextlib
import importlib.util
import io
import sys

from silentgroup import MockGroup, Oracle, load_public
from silentgroup.check import score
from silentgroup.oracle import BUDGET

MOCKS = {
    "pow2": lambda: MockGroup([2**20, 6], seed=0),
    "pow3": lambda: MockGroup([3**18, 3**18 * 7], seed=0),
}

spec = importlib.util.spec_from_file_location("solver", sys.argv[1])
solver = importlib.util.module_from_spec(spec)
spec.loader.exec_module(solver)
solver.VERBOSE = True

name = sys.argv[2] if len(sys.argv) > 2 else "public_1"
runs = int(sys.argv[3]) if len(sys.argv) > 3 else 100
inst = MOCKS[name]() if name in MOCKS else load_public(name)

fails, total = 0, 0.0
for i in range(runs):
    buf = io.StringIO()
    oracle = Oracle(inst, budget=BUDGET)
    with contextlib.redirect_stdout(buf):
        try:
            res = solver.solve(oracle)
        except Exception as e:
            res = None
            print("EXCEPTION", repr(e))
    pts = score(res, inst.answer, BUDGET, oracle.queries_used)
    total += pts
    if pts < 1:
        fails += 1
        print(f"\n===== run {i}: FAILED ({pts:.3f} pts) =====\n{buf.getvalue()}")
    print(f"run {i}: {pts:.3f}", file=sys.stderr)
print(f"\n{fails}/{runs} failed, mean {total / runs:.3f} pts")