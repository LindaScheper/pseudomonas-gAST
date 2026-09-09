"""Abbildung 1: Kohorte, Populationsstruktur und Resistenzlandschaft.

Eingabe : figdata.pkl
Ausgabe : Figure1_cohort.png / .pdf
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

AB = [("CIP","Ciprofloxacin"),("CAZ","Ceftazidime"),("MEM","Meropenem"),("TOB","Tobramycin")]

fig = plt.figure(figsize=(18*CM, 19*CM))
gs = GridSpec(3, 4, figure=fig, height_ratios=[1.45, 1.0, 1.15], hspace=0.6, wspace=0.38)

# --- a) Sequenztypen (gebrochene x-Achse, Top 12 + Other)
# (ehemals Panel b; nimmt jetzt die gesamte obere Reihe ein, da Panel a entfernt wurde)
gs_b = gs[0, :].subgridspec(1, 2, width_ratios=[3.4, 1.0], wspace=0.05)
axL = fig.add_subplot(gs_b[0, 0])
axR = fig.add_subplot(gs_b[0, 1], sharey=axL)

ALL = {s: n for s, n in D["st_all"].items() if s not in ("-", "NA")}
top = sorted(ALL.items(), key=lambda kv: -kv[1])[:12]
n_assigned = sum(ALL.values())
n_other = n_assigned - sum(n for _, n in top)
hr_other = sum(n for s, n in ALL.items()
               if s in D["highrisk"] and s not in dict(top))

rows = [(f"ST{s}", n, s in D["highrisk"]) for s, n in top][::-1]
rows.insert(0, ("Other", n_other, None))

BRK = 300
for ax in (axL, axR):
    for i, (lab, n, hr) in enumerate(rows):
        if lab == "Other":
            ax.barh(i, hr_other, color=OI["vermillion"], height=0.7)
            ax.barh(i, n - hr_other, left=hr_other, color=GREY, height=0.7)
        else:
            ax.barh(i, n, color=OI["vermillion"] if hr else GREY, height=0.7)

axL.set_xlim(0, BRK)
axR.set_xlim(n_other * 0.94, n_other * 1.05)
axR.set_xticks([3100])
axL.set_ylim(-0.7, len(rows) - 0.3)

axL.set_yticks(range(len(rows)))
axL.set_yticklabels([r[0] for r in rows])
axL.tick_params(axis="y", labelsize=7)
axL.tick_params(axis="x", labelsize=7)
axR.tick_params(axis="x", labelsize=7)
axR.tick_params(labelleft=False, left=False)
axR.spines["left"].set_visible(False)
axL.spines["right"].set_visible(False)
axL.set_xlabel("Number of Isolates", fontsize=8); axL.xaxis.set_label_coords(0.68, -0.155)

# Werte an den Balken
for i, (lab, n, _) in enumerate(rows):
    if n > BRK:
        axR.text(n * 1.005, i, f"{n:,}", va="center", ha="left", fontsize=6.5, color="#404040")
    else:
        axL.text(n + 5, i, f"{n:,}", va="center", ha="left", fontsize=6.5, color="#404040")

# Bruchmarkierung auf der Achsenlinie
def _break_marks():
    for ax, at_right in ((axL, True), (axR, False)):
        p = ax.get_position()
        x = p.x1 if at_right else p.x0
        fig.add_artist(plt.Line2D([x - 0.005, x + 0.005], [p.y0 - 0.010, p.y0 + 0.010],
                                  color="k", lw=0.8, clip_on=False))

axL.set_title("Population structure", loc="left")
panel_label(axL, "a", dx=-24)

# --- b) Resistenzlandschaft: Verteilung des Breakpoint-Abstands
# (ehemals Panel c)
gs_c = gs[1, :].subgridspec(1, 5, width_ratios=[1, 1, 1, 1, 0.6], wspace=0.10)
ax_c = [fig.add_subplot(gs_c[0, k]) for k in range(4)]
for k, (ab, name) in enumerate(AB):
    ax = ax_c[k]
    c = D["cohort"][ab]
    bins = np.linspace(-6, 5, 40)
    ax.hist(c["norm_A"], bins=bins, histtype="stepfilled", lw=1.6, color=OI["blue"],
            edgecolor=OI["blue"],
            weights=np.ones(len(c["norm_A"]))/len(c["norm_A"]), label="MIC")
    ax.hist(c["norm_B"], bins=bins, histtype="step", lw=1.3, color=OI["orange"],
            weights=np.ones(len(c["norm_B"]))/len(c["norm_B"]), label="Zone")
    ax.axvline(0, color=OI["vermillion"], lw=1.1, ls="--")
    ax.set_title(name, loc="left")
    ax.set_xlim(-6, 5); ax.set_xticks([-6,-4,-2,0,2,4])
    if k == 0:
        ax.set_ylabel("Fraction of Isolates")
        panel_label(ax, "b", dx=-24)
    ax.tick_params(labelleft=(k == 0))
# Legende rechts neben dem letzten Panel (in den Puffer hinein)
ax_c[3].legend(loc="upper left", bbox_to_anchor=(1.05, 1.0), handlelength=1.1)
ax_c[1].text(1.0, -0.32, "log$_2$ distance from the clinical breakpoint",
             transform=ax_c[1].transAxes, ha="center", fontsize=8)

# --- Reihe 2 (c + d) mit expliziter Pufferspalte dazwischen
# (ehemals d + e)
gs_de = gs[2, :].subgridspec(1, 9, width_ratios=[1,1,1,1, 0.9, 1,1,1,1], wspace=0.15)

x = np.arange(len(AB))
w = 0.30          # etwas schmaler
offset = 0.20     # larger als w/2 -> Gap zwischen den Balken pro Gruppe

# --- c) Resistenzraten
ax = fig.add_subplot(gs_de[0, 0:4])
rA = [D["cohort"][a]["rA"]*100 for a, _ in AB]
rB = [D["cohort"][a]["rB"]*100 for a, _ in AB]
b1 = ax.bar(x-offset, rA, w, color=OI["blue"], label="MIC")
b2 = ax.bar(x+offset, rB, w, color=OI["orange"], label="Zone")
for b in list(b1)+list(b2):
    ax.text(b.get_x()+b.get_width()/2, b.get_height()+0.6, f"{b.get_height():.1f}",
            ha="center", va="bottom", fontsize=6.8)
ax.set_xticks(x); ax.set_xticklabels([n for _, n in AB], rotation=12, ha="right")
ax.set_ylabel("Resistant (%)"); ax.set_ylim(0, 31)
ax.legend(loc="upper right", handlelength=1.2)
ax.set_title("Resistance by antibiotic", loc="left")
panel_label(ax, "c", dx=-24)

# --- d) Isolate vs Patienten
ax = fig.add_subplot(gs_de[0, 5:9])
n_iso = [D["cohort"][a]["n"] for a, _ in AB]
n_pat = [D["npat"][a] for a, _ in AB]
ax.bar(x-offset, n_iso, w, color=GREY, label="Isolates")
ax.bar(x+offset, n_pat, w, color=OI["skyblue"], label="Patients")
for i, (a, b) in enumerate(zip(n_iso, n_pat)):
    ax.text(i-offset, a+150, f"{a:,}", ha="center", va="bottom", fontsize=6.8)
    ax.text(i+offset, b+150, f"{b:,}", ha="center", va="bottom", fontsize=6.8)
ax.set_xticks(x); ax.set_xticklabels([n for _, n in AB], rotation=12, ha="right")
ax.set_ylabel("Count"); ax.set_ylim(0, 8400); ax.set_yticks([0,2000,4000,6000])
ax.legend(loc="upper center", handlelength=1.2, ncol=2, columnspacing=1.4)
ax.set_title("Isolates versus Patients", loc="left")
panel_label(ax, "d", dx=-24)

fig.subplots_adjust(bottom=0.10)
_break_marks()
save(fig, OUT + "SupplementaryFigure1_cohort")