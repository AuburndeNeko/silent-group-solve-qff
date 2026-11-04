# The Silent Group: Track 2 solution

*Qiskit Fall Fest 2026, BITS Pilani. Aditi Kamal*

The challenge itself (problem statement, simulator, checker) is described in `README.md`.
This file describes the solution.

## Quick start

The scored solver is **`my-silent-group/submission.py`** (`solve(oracle) -> dict`; standard
library and qiskit only).

```bash
python -m silentgroup.check my-silent-group/submission.py --repeat 3
```

## Results

Benchmark: both public instances ×10, plus mock groups of six shapes (cyclic and non-cyclic,
~2^20 to ~2^60) × genuine / forged C / forged AB × 5 seeds, scored exactly like the judge with
the 80,000-query budget.

| Version | Correct | Mean score | Mean queries |
|---|---|---|---|
| v2.1 (first robust version) | 81/84 | 1.300 | 25,920 |
| **v4.1 (submitted)** | **84/84** | **1.358** | **22,763** |

Stress tests of v4.1, 100 runs each, all with 0 failures: public_1 (mean 1.346), an
order-2^20 group (1.400) and an order-3^18·7 group (1.382). Slowest run: about 12 s.

## Approach

- **Phase readout:** iterative phase estimation on a single reusable control qubit (measure,
  reset, `if_test` corrections). Each shot samples one character χ, and reads many bits of χ(h).
- **Order of g:** 121-bit readouts, continued fractions, lcm of subsets of readings, device test
  g^N = e. Partial orders backed by repeated readings are lifted by testing small multiples.
- **Exponents:** one joint shot reads χ(g), χ(A), χ(B), χ(C) under the same character. Rounding
  to the 1/n grid gives a·x_g ≡ x_A (mod n). Partial congruences are merged with the CRT, or
  completed by testing their few lifts. Every exponent is verified on the device.
- **Forgery:** a proof via X^n ≠ e when X reads off the grid; otherwise three contradicting shots.
- **Trimming:** deferred. gcd(x_g, N) from the first clean joint shot shows which primes could be
  extra, and is usually 1.

## Open track: forgery forensics

**`open-track/open_track.ipynb`** (all demos run live; start Jupyter from `open-track/`).

A sealed group still leaks its structure, and forging one broadcast doesn't protect the others:

- For each forged broadcast X: the smallest q with X^q = g^w (verified on the device), the
  structure ⟨g, X⟩ ≅ ℤ_d1 × ℤ_d2, a split X = g^c · y into a power of g and a small tag, and how
  several forged tags relate.
- **public_2:** B is forged with tag order q = 585,585; the hidden group is provably
  non-cyclic with |G| ≥ 2.18 × 10¹⁶; the Alice–Charlie pair key is recovered and a message
  encrypted under it is decrypted.
- On mock groups with known structure, the true structure is recovered exactly.
- Also: character-torus visualisations, the solver's circuits next to a photonic
  Mach–Zehnder interferometer, forger fingerprints, case files, and a "be the forger" cell.

The notebook's appendix covers the noise model, the cost floor near 1.36, a lattice (LLL)
experiment that does not pay off under this noise model, and the robustness work.

## Repository layout

```
my-silent-group/
  submission.py      scored solver (v4.1)
  forensics/         open-track code: forensics, decrypt, fingerprints, case files, torus scripts
  testing/           robustness.py (judge-style benchmark), stress.py (repeated runs, logs failures)
  experiments/       lattice self-test and the lattice solver (v3)
  archive/           earlier solver versions, kept for comparison
open-track/          notebook, case files, encrypted message, recorded animation
silentgroup/, examples/   challenge package from the organisers
```

## Reproduce

```bash
python my-silent-group/testing/robustness.py my-silent-group/submission.py --seeds 5 --public 10
python my-silent-group/testing/stress.py my-silent-group/submission.py public_1 100
python my-silent-group/forensics/forensics.py my-silent-group/submission.py
```

The open-track notebook and scripts also need `matplotlib`, `pylatexenc`, `plotly` and
`nbformat`; the submission itself does not.
