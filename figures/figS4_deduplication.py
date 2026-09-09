"""Abbildung 6: Deduplikation nach Patient und Sequenztyp.
Eingabe : figdata.pkl
Ausgabe : Figure6_deduplication.png / .pdf
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
ABS  = ["CIP", "MEM", "CAZ", "TOB"]   # Reihenfolge wie in der Abbildung
D = pickle.load(open(DATA + "figdata.pkl", "rb"))
# =========================== Abbildung 6 ===========================
fig = plt.figure(figsize=(18*CM, 8.5*CM))
gs = GridSpec(1, 3, figure=fig, wspace=0.44, width_ratios=[1, 1, 1.15])
# a) Praevalenz und BSS unter Deduplikation
ax = fig.add_subplot(gs[0, 0])
x = np.arange(len(ABS)); w = 0.26
lab3 = ["full collection", "first isolate", "most resistant"]
sh = [GREY, OI["skyblue"], OI["orange"]]
bars_for_legend = []
for k in range(3):
    b_ = ax.bar(x + (k-1)*w, [D["dedup"][a]["prev"][k] for a in ABS], w, color=sh[k], label=lab3[k])
    bars_for_legend.append(b_)
ax.set_xticks(x); ax.set_xticklabels([NAME[a] for a in ABS], rotation=22, ha="right")
ax.set_ylabel("Resistant (%)"); ax.set_ylim(0, 24)
ax.set_title("Prevalence", loc="left", fontsize=8.5)
panel_label(ax, "a", dx=-24.0)
ax_a = ax
# b) Brier Skill Score
ax = fig.add_subplot(gs[0, 1])
for k in range(3):
    ax.bar(x + (k-1)*w, [D["dedup"][a]["bss"][k] for a in ABS], w, color=sh[k])
ax.axhline(0, color="k", lw=0.7)
ax.set_xticks(x); ax.set_xticklabels([NAME[a] for a in ABS], rotation=22, ha="right")
ax.set_ylabel("Brier skill score"); ax.set_ylim(0, 0.62)
ax.set_title("Predictive skill", loc="left", fontsize=8.5)
panel_label(ax, "b", dx=-24.0)
ax_b = ax
# c) Welche Marker ueberleben
ax = fig.add_subplot(gs[0, 2])
# nachgerechnet in vorstufen/13_dedup_marker.py; die frueher hier verwendete
# feste Liste in 01_figdata.py war fuer PBP4 und fusA1 nicht mehr aktuell
dm = [(r[0], r[1], r[2], r[3]) for r in json.load(open(DATA + "dedup_markers.json"))][::-1]
for i, (m, full, dd, n) in enumerate(dm):
    c = OI["green"] if dd/full >= 0.80 else OI["vermillion"]   # >=80 % des urspruenglichen PPV
    ax.plot([full, dd], [i, i], color=LGREY, lw=1.2, zorder=1)
    ax.plot(full, i, "o", color=GREY, ms=4, zorder=2)
    ax.plot(dd, i, "o", color=c, ms=5.5, zorder=3)
ax.set_yticks(range(len(dm)))
def _fmt_label(m):
    parts = m.split(maxsplit=1)
    if len(parts) == 2:
        gene, rest = parts
        return f"$\\it{{{gene}}}$ {rest}"
    return f"$\\it{{{m}}}$"

ax.set_yticklabels([_fmt_label(m) for m, _, _, _ in dm], fontsize=6.6)
ax.set_xlabel("Positive predictive value (%)"); ax.set_xlim(20, 105)
ax.plot([], [], "o", color=GREY, ms=4, label="full collection")
ax.plot([], [], "o", color=OI["green"], ms=5, label="de-duplicated, retained")
ax.plot([], [], "o", color=OI["vermillion"], ms=5, label="de-duplicated, lost")
ax.legend(fontsize=6.3, handlelength=1.0, loc="upper left",
          bbox_to_anchor=(-0.02, -0.20), ncol=1)
ax.set_title("Marker retention", loc="left", fontsize=8.5)
panel_label(ax, "c", dx=-24.0)
ax_c = ax

# ---- Abstaende a-b und b-c separat justieren ------------------------------
# gap_ab / gap_bc: positiver Wert = groesserer Abstand, negativer Wert = kleinerer Abstand
# (Werte sind Anteile der Figure-Breite, typischer Bereich etwa -0.03 .. +0.05)
gap_ab = 0.0
gap_bc = 0.05

pos_a = ax_a.get_position()
pos_b = ax_b.get_position()
pos_c = ax_c.get_position()

shift_b = gap_ab
shift_c = gap_ab + gap_bc

ax_b.set_position([pos_b.x0 + shift_b, pos_b.y0, pos_b.width, pos_b.height])
ax_c.set_position([pos_c.x0 + shift_c, pos_c.y0, pos_c.width, pos_c.height])

# gemeinsame Legende fuer a und b, mittig darunter, nebeneinander
# (Positionen nach der Verschiebung neu abfragen!)
pos_a = ax_a.get_position(); pos_b = ax_b.get_position()
center_x = (pos_a.x0 + pos_b.x1) / 2
bottom_y = min(pos_a.y0, pos_b.y0) - 0.14
fig.legend(bars_for_legend, lab3, loc="upper center",
           bbox_to_anchor=(center_x, bottom_y), bbox_transform=fig.transFigure,
           ncol=3, fontsize=6.5, handlelength=1.0, frameon=False)

save(fig, OUT + "SupplementaryFigure4_deduplication")