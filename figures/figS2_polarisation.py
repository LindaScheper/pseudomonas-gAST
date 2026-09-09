"""Supplementary Figure 2: Referenzpolarisierung und Pseudomarker.
Eingabe : figdata.pkl
Ausgabe : SupplementaryFigure2_polarisation.png / .pdf
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
fig = plt.figure(figsize=(18*CM, 12*CM))
gs = GridSpec(2, 3, figure=fig, hspace=0.72, wspace=0.40, width_ratios=[1, 1, 0.005])
# a) Frequenzspektrum
ax = fig.add_subplot(gs[0, 0])
v = D["varfreq"]
ax.hist(v, bins=np.linspace(0, 1, 51), color=GREY, edgecolor="white", linewidth=0.3)
ax.axvspan(0.90, 1.0, color=OI["vermillion"], alpha=0.16)
ax.set_yscale("log")
ax.set_xlim(0, 1)
ax.set_xlabel("Population frequency of the variant allele")
ax.set_ylabel("Variant positions (n)")
ax.set_title("Variant frequency spectrum", loc="left", fontsize=8.5)
panel_label(ax, "a", dx=-24.0)
# b) Die haeufigsten Pseudo-Marker
ax = fig.add_subplot(gs[0, 1:])
nf = D["nearfixed"][:9]
items = [(g, m, f*100) for g, p, m, f in nf] + [("nalC", "G71E", 91.3)]
items = items[::-1]
lbl = [f"$\\it{{{g}}}$ {m}" for g, m, _ in items]
val = [v_ for _, _, v_ in items]
for i_, v_ in enumerate(val):
    ax.plot([88, v_], [i_, i_], color=LGREY, lw=0.9, zorder=1)
    ax.plot(v_, i_, "o", color=GREY, ms=6, zorder=2)
    ax.text(v_ + 0.35, i_, f"{v_:.1f}", va="center", fontsize=6.5, color="black")
ax.set_yticks(range(len(items))); ax.set_yticklabels(lbl, fontsize=7)
ax.set_xlim(88, 102); ax.set_xlabel("Isolates carrying the allele (%)")
ax.set_title("Near-fixed variants misclassified as markers", loc="left", fontsize=8.5)
pos = ax.get_position()
gap = 0.05
ax.set_position([pos.x0 + gap, pos.y0, pos.width - gap, pos.height])
panel_label(ax, "b", dx=-24.0)
# c) Spanische Gegenprobe
ax = fig.add_subplot(gs[1, :2])
sp = D["spain_filter"]
x = np.arange(len(sp))
rem = [(a-b)/a*100 for _, a, b in sp]
ax.bar(x, rem, color=[OI["skyblue"] if r < 60 else OI["orange"] for r in rem], width=0.68)
for i_, r in enumerate(rem):
    ax.text(i_, r+1.5, f"{r:.0f}", ha="center", fontsize=6.5)
ax.set_xticks(x); ax.set_xticklabels([g for g, _, _ in sp], rotation=45, ha="right", style="italic")
ax.set_ylabel("Variantallels removed (%)"); ax.set_ylim(0, 108)
ax.set_title("Validation in an independent cohort", loc="left", fontsize=8.5)
panel_label(ax, "c", dx=-24.0)
save(fig, OUT + "SupplementaryFigure2_polarisation")