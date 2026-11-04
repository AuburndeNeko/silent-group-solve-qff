"""Case files: one dossier per meeting, generated from fresh forensics runs.

Saves casefiles.html (open it in a browser). Includes fingerprint.png and torus3d.png
if they exist in the current folder.

    python my-silent-group/casefiles.py my-silent-group/submission.py
"""

import datetime
import html
import importlib.util
import math
import os
import sys

here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("forensics", os.path.join(here, "forensics.py"))
F = importlib.util.module_from_spec(spec)
spec.loader.exec_module(F)          # uses sys.argv[1] as the solver, like forensics.py

CSS = """
body { background:#e9e4d8; font-family:'Courier New',monospace; color:#222; margin:0; padding:30px; }
.case { background:#fbf8f0; max-width:860px; margin:0 auto 40px; padding:30px 40px;
        border:1px solid #b9b2a0; box-shadow:4px 4px 0 #cfc8b6; position:relative; }
h1 { font-size:22px; letter-spacing:2px; border-bottom:2px solid #222; padding-bottom:6px; }
h2 { font-size:15px; letter-spacing:1px; margin-top:22px; text-transform:uppercase; }
.stamp { position:absolute; top:28px; right:36px; transform:rotate(-8deg); font-size:20px;
         font-weight:bold; padding:6px 14px; border:4px solid; letter-spacing:2px; }
.forged { color:#b3261e; border-color:#b3261e; }
.auth { color:#1e7b34; border-color:#1e7b34; }
table { width:100%; border-collapse:collapse; font-size:13px; }
td, th { border-bottom:1px dashed #aaa; padding:6px; vertical-align:top; text-align:left; }
.tag { font-weight:bold; }
.g { color:#1e7b34; } .r { color:#b3261e; }
.note { font-size:12px; color:#666; }
img { max-width:100%; border:1px solid #ccc; margin-top:8px; }
"""


def esc(s):
    return html.escape(str(s))


def case_html(num, label, summary):
    n = summary["n"]
    forged = [X for X in "ABC" if summary[X][0] == "forged"]
    genuine = {X: summary[X][1] for X in "ABC" if summary[X][0] == "genuine"}
    stamp = ('<div class="stamp forged">FORGERY DETECTED</div>' if forged
             else '<div class="stamp auth">AUTHENTIC</div>')
    rows = []
    for X in "ABC":
        info = summary[X]
        if info[0] == "genuine":
            rows.append(f"<tr><td class='tag'>{X}</td><td class='g'>AUTHENTIC</td>"
                        f"<td>{X} = g<sup>{info[1]}</sup><br><span class='note'>verified on the device</span></td></tr>")
        elif len(info) == 3:
            rows.append(f"<tr><td class='tag'>{X}</td><td class='r'>FORGED</td>"
                        f"<td>ord({X}) = {info[1]}<br><span class='note'>relation to &lang;g&rang; not determined</span></td></tr>")
        else:
            _, nX, q, w, d1, d2, split = info
            ev = [f"ord({X}) = {nX}",
                  f"smallest q with {X}<sup>q</sup> &isin; &lang;g&rang;: <b>q = {q}</b>",
                  f"verified on device: {X}<sup>{q}</sup> = g<sup>{w}</sup>",
                  f"&lang;g, {X}&rang; &cong; Z<sub>{d1}</sub> &times; Z<sub>{d2}</sub>"]
            if split:
                ev.append(f"split: {X} = g<sup>{split[0]}</sup> &middot; y, with y<sup>{q}</sup> = e "
                          f"<span class='note'>(c defined mod {split[1]})</span>")
            rows.append(f"<tr><td class='tag'>{X}</td><td class='r'>FORGED</td><td>{'<br>'.join(ev)}</td></tr>")

    keys = [(k[4:], v) for k, v in summary.items() if k.startswith("key_")]
    key_lines = [f"pair key g<sup>({p.lower()})</sup> = g<sup>{v}</sup> "
                 f"<span class='note'>(checked on the device)</span>" for p, v in keys]
    if len(genuine) == 3:
        key_lines.append(f"<b>full three-party secret abc mod n = "
                         f"{genuine['A'] * genuine['B'] * genuine['C'] % n}</b>")
    if not key_lines:
        key_lines = ["no pair of genuine parties: no keys recoverable"]

    facts = []
    if "Q" in summary:
        order_Q, expQ, rels = summary["Q"]
        facts.append(f"tag group Q = &lang;g, A, B, C&rang; / &lang;g&rang;: |Q| = {order_Q}"
                     + ("" if expQ == order_Q else f", exponent {expQ}"))
        facts.append("tags " + ("share relations: the same fake was reused" if rels
                                else "are independent: different fakes"))
        bound = n * order_Q
    elif len(forged) == 1 and len(summary[forged[0]]) > 3:
        bound = n * summary[forged[0]][2]
    else:
        bound = n
    facts.append(f"|G| &ge; {bound} (&asymp; 2<sup>{math.log2(bound):.1f}</sup>)")
    nc = [X for X in forged if len(summary[X]) > 3 and summary[X][4] > 1]
    if nc:
        X = nc[0]
        d1, d2 = summary[X][4], summary[X][5]
        facts.append(f"<b>G is not cyclic</b>: &lang;g, {X}&rang; &cong; Z<sub>{d1}</sub> &times; "
                     f"Z<sub>{d2}</sub> with gcd {d1} &gt; 1, and every subgroup of a cyclic group is cyclic")

    verdict = (f"Broadcast(s) {', '.join(forged)} forged. No three-party key exists"
               + (f"; the honest pair's channel remains readable." if keys else ".")
               if forged else "All three broadcasts authentic. Three-party key recovered.")
    return f"""
<div class="case">{stamp}
<h1>CASE FILE #{num:03d}: {esc(label)}</h1>
<p><b>Order of the public generator:</b> n = ord(g) = {n}</p>
<h2>Parties</h2>
<table><tr><th>Broadcast</th><th>Verdict</th><th>Evidence</th></tr>{''.join(rows)}</table>
<h2>Recovered keys</h2><p>{'<br>'.join(key_lines)}</p>
<h2>About the hidden group</h2><p>{'<br>'.join(facts)}</p>
<h2>Verdict</h2><p>{verdict}</p>
</div>"""

if __name__ == "__main__":
    cases = [c for c in F.CASES if not c[0].startswith("mock non-cyclic")] + \
            [c for c in F.CASES if c[0].startswith("mock non-cyclic")]
    pages = []
    for i, (label, make, moduli) in enumerate(cases, start=1):
        print(f"investigating {label} ...", flush=True)
        pages.append(case_html(i, label, F.investigate(label, make(), moduli, quiet=True)))

    extras = ""
    for img, cap in (("fingerprint.png", "Forger fingerprints: prime factors of each tag order"),
                     ("torus3d.png", "Character torus (illustration on a small mock group)")):
        if os.path.exists(img):
            extras += f"<h2>{cap}</h2><img src='{img}'>"
    if extras:
        pages.append(f"<div class='case'><h1>EXHIBITS</h1>{extras}</div>")

    with open("casefiles.html", "w") as f:
        f.write(f"<!doctype html><html><head><meta charset='utf-8'><title>Case files</title>"
                f"<style>{CSS}</style></head><body>"
                f"<p class='note' style='text-align:center'>The Silent Group: forensic case files, "
                f"generated {datetime.date.today()}</p>{''.join(pages)}</body></html>")
    print("saved casefiles.html")