"""Partial key recovery, made concrete: read the honest parties' private channel.

In public_2, B's broadcast is forged, so there is no three-party key. But A and C are
genuine, so they still share the pair key g^(ac) and can use it to talk privately.

  1) Honest side (run once): encrypt a message under the pair-key exponent ac mod n.
       python my-silent-group/decrypt.py --make "your message" 23054424601
  2) Attacker side: recover a and c on the device from scratch, then decrypt.
       python my-silent-group/decrypt.py my-silent-group/submission.py
"""

import hashlib
import hmac
import importlib.util
import json
import os
import sys

FILE = "message.enc"
PAIR = "AC"


def derive_keys(k):
    base = hashlib.sha256(b"silent-group pair key:" + str(k).encode()).digest()
    return hashlib.sha256(base + b"enc").digest(), hashlib.sha256(base + b"mac").digest()


def keystream(key, nonce, length):
    out, i = b"", 0
    while len(out) < length:
        out += hashlib.sha256(key + nonce + i.to_bytes(8, "big")).digest()
        i += 1
    return out[:length]


def encrypt(msg, k):
    ek, mk = derive_keys(k)
    nonce = os.urandom(16)
    data = msg.encode()
    ct = bytes(x ^ y for x, y in zip(data, keystream(ek, nonce, len(data))))
    tag = hmac.new(mk, nonce + ct, hashlib.sha256).hexdigest()
    return {"pair": PAIR, "nonce": nonce.hex(), "ct": ct.hex(), "tag": tag}


def decrypt(box, k):
    ek, mk = derive_keys(k)
    nonce, ct = bytes.fromhex(box["nonce"]), bytes.fromhex(box["ct"])
    expected = hmac.new(mk, nonce + ct, hashlib.sha256).hexdigest()
    if not hmac.compare_digest(expected, box["tag"]):
        return None
    return bytes(x ^ y for x, y in zip(ct, keystream(ek, nonce, len(ct)))).decode()


if sys.argv[1] == "--make":
    with open(FILE, "w") as f:
        json.dump(encrypt(sys.argv[2], int(sys.argv[3])), f, indent=1)
    print(f"wrote {FILE} (ciphertext + tag only; the key is not stored)")
    sys.exit()

from silentgroup import Oracle, load_public

spec = importlib.util.spec_from_file_location("solver", sys.argv[1])
S = importlib.util.module_from_spec(spec)
spec.loader.exec_module(S)
S.VERBOSE = False

with open(FILE) as f:
    box = json.load(f)
print(f"intercepted ciphertext: {box['ct'][:48]}...")
print("recovering the broadcasts of public_2 from scratch on the device:")

oracle = Oracle(load_public("public_2"))
n = S.find_order(oracle, "g", trim=True)
exps = {}
for X in S.ELEMENTS:
    found, _, _ = S.find_exponents(oracle, n, elements=(X,))
    if X in found:
        exps[X] = found[X]
        print(f"  {X} = g^{found[X]}   (genuine, verified)")
    else:
        print(f"  {X}: forged -> no three-party key exists")

X, Y = box["pair"]
if X in exps and Y in exps:
    k = exps[X] * exps[Y] % n
    print(f"pair key exponent {X.lower()}{Y.lower()} mod n = {k}")
    msg = decrypt(box, k)
    print(f"decrypted: {msg!r}" if msg is not None else "decryption FAILED (tag mismatch)")
    wrong = decrypt(box, k + 1)
    print(f"same ciphertext, wrong key: {'rejected (tag mismatch)' if wrong is None else wrong!r}")
else:
    print(f"pair {X}{Y} is not fully genuine: cannot open this channel")
print(f"(queries used: {oracle.queries_used})")