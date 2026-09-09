"""Abbildung 2: Antibiotikaspezifitaet auf der Odds-Ratio-Skala, eine Tafel.

Elf benannte Bloecke mit ihrem Spezifitaetsquotienten und Konfidenzintervall auf
der Odds-Ratio-Skala (Punkt) und derselbe Block auf der Vorhersagewert-Skala
(graues Kreuz). Farbig, wenn der Block das Kriterium erfuellt, sonst grau.

Das Nullexperiment, das frueher Tafel a war, steht jetzt in
figS3_nullexperiment.py; die Streuwolke der Umklassifikation ist entfallen, ihre
beiden Zahlen stehen im Text.

Eingabe : specratio.pkl, specratio.json
Ausgabe : Figure2_specificity.png / .pdf
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import *
from paths import DATA, OUT
import json
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

COL  = {"CIP": OI["vermillion"], "CAZ": OI["green"], "MEM": OI["blue"], "TOB": OI["purple"]}
ABS  = ["CIP", "CAZ", "MEM", "TOB"]

R = pd.read_pickle(DATA + "specratio.pkl")

sel = [("gyrA_83", "CIP", "$\\it{gyrA}$ T83 → CIP"),
       ("parC_87", "CIP", "$\\it{parC}$ S87 → CIP"),
       ("TRUNC:PBP4", "CAZ", "$\\it{PBP4}$ truncated → CAZ"),
       ("TRUNC:ampD", "CAZ", "$\\it{ampD}$ truncated → CAZ"),
       ("fusA1_671", "TOB", "$\\it{fusA1}$ A671 → TOB"),
       ("TRUNC:mexA", "CIP", "$\\it{mexA}$ truncated → CIP"),
       ("TRUNC:oprD", "MEM", "$\\it{oprD}$ truncated → MEM"),
       ("LOSS:oprD", "MEM", "$\\it{oprD}$ absent → MEM"),
       ("ELONG:oprD", "MEM", "$\\it{oprD}$ elongated → MEM"),
       ("STOP:mexR", "MEM", "$\\it{mexR}$ nonsense → MEM"),
       ("LOSS:mexZ", "TOB", "$\\it{mexZ}$ absent → TOB"),
       ("mexZ_126", "TOB", "$\\it{mexZ}$ clone block → TOB"),
       ("ampO_45 [+62]", "CIP", "$\\it{ampO}$ block → CIP")]
rows = []
for mk, ab, lab in sel:
    q = R[(R.marker == mk) & (R.ab == ab)]
    if len(q):
        r = q.iloc[0]
        rows.append((lab, ab, r.ratio, r.lo, r.hi, r.ratio_ppv, int(r.n), bool(r.spec)))
rows.sort(key=lambda t: t[2])

fig, ax = plt.subplots(figsize=(13*CM, 9.5*CM))
for i, (lab, ab, v, lo, hi, rp, n, sp) in enumerate(rows):
    c = COL[ab] if sp else GREY
    ax.plot([lo, hi], [i, i], color=c, lw=1.5, solid_capstyle="round", alpha=0.85, zorder=2)
    ax.plot(v, i, "o", color=c, ms=6, zorder=3)
    ax.plot(rp, i, "x", color="#B0B0B0", ms=5, mew=1.3, zorder=1)
    ax.text(1.02, i, f"n = {n}", transform=ax.get_yaxis_transform(),
            va="center", fontsize=6.4, color="#606060")
ax.axvline(1.4, color="k", lw=1.0, ls="--")
ax.axvline(1.0, color=LGREY, lw=0.8, zorder=0)
ax.set_xscale("log")
ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[0] for r in rows], fontsize=7.4)
ax.set_xlim(0.5, 400)
ax.set_xticks([1, 1.4, 10, 100])
ax.set_xticklabels(["1", "1.4", "10", "100"])
ax.set_xlabel("Antibiotic-specificity ratio (95 % CI)\n"
              "Odds Ratio for this antibiotic/geometric mean for the other three")
# Legende für Punkt/Kreuz
legend_scale = [
    Line2D([0], [0], marker="o", color=GREY,
           markersize=5, linestyle="none",
           label="Odds ratio scale"),
    Line2D([0], [0], marker="x", color="#B0B0B0",
           markersize=5, markeredgewidth=1.3,
           linestyle="none",
           label="predictive-value scale")
]

# Legende für Antibiotika-Farben
legend_ab = [
    Line2D([0], [0], color=COL["CIP"], marker="o", linestyle="-",
           markersize=5, linewidth=1.2, label="Ciprofloxacin"),
    Line2D([0], [0], color=COL["CAZ"], marker="o", linestyle="-",
           markersize=5, linewidth=1.2, label="Ceftazidime"),
    Line2D([0], [0], color=COL["MEM"], marker="o", linestyle="-",
           markersize=5, linewidth=1.2, label="Meropenem"),
    Line2D([0], [0], color=COL["TOB"], marker="o", linestyle="-",
           markersize=5, linewidth=1.2, label="Tobramycin")
]

leg1 = ax.legend(
    handles=legend_scale,
    loc="upper center",
    bbox_to_anchor=(0.5, -0.22),
    fontsize=7.5,
    ncol=2.0,
    handlelength=1.4,
    columnspacing=1.2,
    labelspacing=0.3
)

leg2 = ax.legend(
    handles=legend_ab,
    loc="upper center",
    bbox_to_anchor=(0.5, -0.31),
    fontsize=7.5,
    ncol=4.0,
    handlelength=1.4,
    columnspacing=1.2,
    labelspacing=0.3
)

ax.add_artist(leg1)
fig.subplots_adjust(left=0.30, right=0.88)
save(fig, OUT + "Figure2_specificity")
