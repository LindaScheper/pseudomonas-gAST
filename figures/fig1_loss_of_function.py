"""Abbildung 2: Proteinlaenge und explizite Inaktivierungsaufrufe.

Eingabe : figdata.pkl, fs_fig.json
Ausgabe : Figure1_loss_of_function.png / .pdf
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import *
from paths import DATA, OUT
import pickle, json, collections
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

COL  = {"CIP": OI["vermillion"], "CAZ": OI["green"], "MEM": OI["blue"], "TOB": OI["purple"]}
NAME = {"CIP": "Ciprofloxacin", "CAZ": "Ceftazidime", "MEM": "Meropenem", "TOB": "Tobramycin"}
ABS  = ["CIP", "CAZ", "MEM", "TOB"]

D = pickle.load(open(DATA + "figdata.pkl", "rb"))

by_m = collections.defaultdict(dict)
for t in D["trunc"]:
    by_m[t["m"]][t["ab"]] = t
rows = [(m, v) for m, v in by_m.items() if max(t["n"] for t in v.values()) >= 20]
rows.sort(key=lambda r: -max(t["orr"] for t in r[1].values()))
# Die zehn staerksten Assoziationen insgesamt, ergaenzt um die zwei staerksten
# jeder Substanz, die darunter sonst nicht vorkaeme. Ohne diese Ergaenzung
# fehlte Tobramycin ganz: sein staerkster Verlustmarker (mexZ) erreicht nur
# eine Odds Ratio von 4,4 und faellt hinter die Beta-Laktam-Loci zurueck.
keep = [m for m, _ in rows[:10]]
for ab in ABS:
    if any(ab in by_m[m] for m in keep):
        continue
    c = sorted([(m, v[ab]) for m, v in by_m.items() if ab in v and v[ab]["n"] >= 20],
               key=lambda x: -x[1]["orr"])
    keep += [m for m, _ in c[:2] if m not in keep]
rows = [(m, by_m[m]) for m in keep]
rows.sort(key=lambda r: -max(t["orr"] for t in r[1].values()))
rows = rows[::-1]

fig, ax = plt.subplots(figsize=(13*CM, 8.5*CM))
for i_, (m, v) in enumerate(rows):
    ax.axhline(i_, color=LGREY, lw=0.6, zorder=0)
    for ab in ["CIP", "CAZ", "MEM", "TOB"]:
        if ab not in v: continue
        t = v[ab]
        ax.plot([t["lo"], t["hi"]], [i_, i_], color=COL[ab], lw=1.2, alpha=0.8,
                solid_capstyle="round", zorder=2)
        ax.plot(t["orr"], i_, "o", color=COL[ab], ms=5.5, zorder=3)
ax.axvline(1, color="k", lw=0.9, ls=":")
ax.set_xscale("log"); ax.set_yticks(range(len(rows)))
def nice(m):
    g = m.split(":")[1]
    return f"$\\it{{{g}}}$ {'truncated' if m.startswith('TRUNC') else 'absent'}"
ax.set_yticklabels([nice(m) for m, _ in rows], fontsize=8)
ax.set_xlabel("Odds Ratio for resistance,\ncarriers versus non-carriers (95 % CI)")
ax.set_xlim(0.8, 200)
# Traegerzahl und Resistenzrate der staerksten Substanz rechts anschreiben,
# in der Farbe des jeweiligen Antibiotikums
for i_, (m, v) in enumerate(rows):
    ab = max(v, key=lambda a: v[a]["orr"]); t = v[ab]
    ax.text(1.03, i_, f"n = {t['n']},  {t['ppv']*100:.0f} %",
            transform=ax.get_yaxis_transform(), va="center", fontsize=6.6, color=COL[ab])
# x-Grenze dynamisch: alle OR-Punkte > 1 der gezeigten Zeilen muessen sichtbar sein
max_orr = max(t["orr"] for _, v in rows for t in v.values())
ax.set_xlim(0.6, max_orr * 1.15)
for ab in ["CIP", "CAZ", "MEM", "TOB"]:
    if any(ab in v for _, v in rows): ax.plot([], [], "o-", color=COL[ab], label=NAME[ab], ms=4.5)
ax.legend(loc="upper center", fontsize=7.5, handlelength=1.4, ncol=4.0,
          bbox_to_anchor=(0.5, -0.22), columnspacing=1.2)
fig.subplots_adjust(right=0.62)
save(fig, OUT + "Figure1_loss_of_function")