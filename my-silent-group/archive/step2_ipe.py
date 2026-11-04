"""Step 2: read chi(g) to m bits in ONE shot, using one reusable control qubit."""

import math
from fractions import Fraction

from qiskit import ClassicalRegister, QuantumCircuit, QuantumRegister

from silentgroup import MockGroup, Oracle, oracle_gate


def ipe_circuit(element: str, m: int) -> QuantumCircuit:
    """Iterative phase estimation of chi(element) to m bits.

    Bit b_k (k = 1..m, with phi = 0.b1 b2 ... bm) is stored in classical bit k-1.
    """
    q = QuantumRegister(1, "q")
    c = ClassicalRegister(m, "c")
    qc = QuantumCircuit(q, c)
    for k in range(m, 0, -1):              # biggest power first: k = m, m-1, ..., 1
        if k != m:
            qc.reset(q[0])                 # reuse the same qubit
        qc.h(q[0])
        qc.append(oracle_gate(element, 2 ** (k - 1)), [q[0]])   # rotation = 0.b_k b_(k+1)... turns
        for l in range(k + 1, m + 1):      # undo the bits we already know
            with qc.if_test((c[l - 1], 1)):
                qc.p(-math.pi / 2 ** (l - k), q[0])
        qc.h(q[0])                         # now the rotation is 0 or 1/2 turn -> measure reads b_k
        qc.measure(q[0], c[k - 1])
    return qc


def counts_to_phases(counts: dict, m: int) -> list:
    """Exact phase estimates as Fractions (never floats: floats keep only ~53 bits)."""
    phases = []
    for key, n in counts.items():
        phases += [Fraction(int(key[::-1], 2), 2 ** m)] * n
    return phases


def experiment(label, mock, m, shots, noise, max_den):
    oracle = Oracle(mock, noise=noise) if noise is not None else Oracle(mock)
    qc = ipe_circuit("g", m)
    print(f"\n=== {label}: m={m} bits, cost per shot = {oracle.cost(qc, 1)} ===")
    print("answer:", mock.answer)
    counts = oracle.run(qc, shots=shots)
    for phi in counts_to_phases(counts, m):
        f = Fraction(phi).limit_denominator(max_den)
        print(f"  phi ~ {phi:.12f}   best fraction s/r = {f}   r = {f.denominator}")
    print("queries used:", oracle.queries_used)


if __name__ == "__main__":
    # A: Z_256, phases are exact multiples of 1/256, so 8 bits should be exact.
    experiment("A (exact, no noise)", MockGroup([256], g=[3], seed=0),
               m=8, shots=10, noise=0, max_den=256)

    # B: Z_1000003 (prime ~2^20). Phases are s/1000003: not exact in binary.
    #    41 bits ~ 2 x 20 bits, enough for continued fractions to recover r.
    mock_b = MockGroup([1000003], g=[5], seed=0)
    experiment("B (inexact, no noise)", mock_b, m=41, shots=10, noise=0, max_den=2**20)

    # C: same as B, with the device's default (realistic) noise.
    experiment("C (inexact, default noise)", mock_b, m=41, shots=10, noise=None, max_den=2**20)