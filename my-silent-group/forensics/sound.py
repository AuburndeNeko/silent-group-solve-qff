"""Forgeries sound out of tune.

Each shot -> a two-note chord: the PREDICTED phase a*chi(g) and the ACTUAL phase chi(X),
as pitches in one octave (phase p -> 330 * 2^p Hz). Genuine X = g^a: always a clean
unison. Forged X: no a works, so the notes clash. Writes genuine.wav and forged.wav.

    python my-silent-group/sound.py my-silent-group/submission.py
"""

import importlib.util
import sys
import wave

import numpy as np

from silentgroup import MockGroup, Oracle

spec = importlib.util.spec_from_file_location("solver", sys.argv[1])
S = importlib.util.module_from_spec(spec)
spec.loader.exec_module(S)
S.VERBOSE = False

M, SHOTS, BITS = 31, 24, 14
RATE, BASE, NOTE = 44100, 330.0, 0.3      # sample rate, lowest pitch (Hz), seconds per chord

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
    shots.append({k: round(v * M) % M for k, v in ph.items()})    # grid points x/M


def guess_exponent(X):
    """Exponent read off the first shot with chi(g) != 0 (exact if X is genuine)."""
    for s in shots:
        if s["g"]:
            return s[X] * pow(s["g"], -1, M) % M
    return 0


def tone(f, n):
    t = np.arange(n) / RATE
    env = np.minimum(1.0, np.minimum(t, t[::-1]) / 0.02)           # 20 ms fade in/out
    return np.sin(2 * np.pi * f * t) * env


def render(X, fname):
    a = guess_exponent(X)
    n = int(RATE * NOTE)
    chunks, cents = [], []
    for s in shots:
        predicted = (a * s["g"] % M) / M
        actual = s[X] / M
        f1, f2 = BASE * 2 ** predicted, BASE * 2 ** actual
        cents.append(abs(1200 * np.log2(f2 / f1)))
        chunks.append(0.4 * (tone(f1, n) + tone(f2, n)))
    audio = np.concatenate(chunks)
    with wave.open(fname, "wb") as w:
        w.setnchannels(1)
        w.setsampwidth(2)
        w.setframerate(RATE)
        w.writeframes((audio * 32767).astype(np.int16).tobytes())
    in_tune = sum(c < 1 for c in cents)
    print(f"{fname}: {X} with guessed exponent {a}: {in_tune}/{SHOTS} chords in tune, "
          f"mean detune {np.mean(cents):.0f} cents")


render("A", "genuine.wav")
render("C", "forged.wav")