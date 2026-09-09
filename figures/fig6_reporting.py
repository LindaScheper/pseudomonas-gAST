"""Figure 6: Befundzustaende und Kalibrierung.

Tafel a  Anteil der Isolate je Befundzustand (nicht-resistent / blank /
         resistent), oben aus der Modellwahrscheinlichkeit allein, unten mit
         dem Enzym-Rule-in.
Tafel b  Beobachtete Resistenz in den beiden berichteten Gruppen, gegen die
         VME-Schranke von 3 % und die Rule-in-Schranke von 80 %.
Tafel c  Lage der drei Befundzustaende auf der MIC-Skala (Zentrum A), als
         Median und Interquartilsabstand des log2-Abstands vom Grenzwert.

Die Zahlen fuer a und b kommen aus rulein_v2.json, dem Lauf mit dem
vorzeichenbehafteten Pumpenmerkmal, aus dem auch die Zahlen im Manuskript
stammen. Das aeltere D["rulein"] in figdata.pkl wird nicht mehr verwendet.
Tafel c kommt aus drei_zustaende.json.

Eingabe : rulein_v2.json, drei_zustaende.json
Ausgabe : Figure6_reporting.png / .pdf
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import *
from paths import DATA, OUT
import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch

COL  = {"CIP": OI["vermillion"], "CAZ": OI["green"], "MEM": OI["blue"], "TOB": OI["purple"]}
NAME = {"CIP": "Ciprofloxacin", "CAZ": "Ceftazidime", "MEM": "Meropenem", "TOB": "Tobramycin"}
ABS  = ["CIP", "MEM", "CAZ", "TOB"]   # Reihenfolge wie in der Tabelle im Manuskript

RI = json.load(open(DATA + "rulein_v2.json"))
# Die Datei haelt Anteile; die Abbildung rechnet in Prozent.
PCT = ["S", "B", "Sr", "R_old", "R_new", "Rr_old", "Rr_new"]
for ab in ABS:
    for k in PCT:
        if RI[ab].get(k) is not None:
            RI[ab][k] = RI[ab][k] * 100

fig = plt.figure(figsize=(21*CM, 8.5*CM))
gs = GridSpec(1, 3, figure=fig, wspace=0.48, width_ratios=[1.15, 1, 1.0])

# --- a) Befundzustaende, mit und ohne Enzym-Rule-in -------------------------
ax = fig.add_subplot(gs[0, 0])
y = np.arange(len(ABS))[::-1]
h = 0.34
MIN_W = 8.0  # unterhalb dieser Segmentbreite (%) wird die Zahl nach aussen verlegt
for i, ab in enumerate(ABS):
    v = RI[ab]
    S = v["S"]
    for key, off, alpha in [("R_old", +h/2, 0.45), ("R_new", -h/2, 1.0)]:
        R_ = v[key]
        B  = 100 - S - R_
        ax.barh(y[i]+off, S,  color=OI["skyblue"],    height=h, alpha=alpha)
        ax.barh(y[i]+off, B,  left=S, color=LGREY,    height=h, alpha=alpha)
        ax.barh(y[i]+off, R_, left=S+B, color=OI["vermillion"], height=h, alpha=alpha)

        # S liegt links im Balken -> unter 1% keine Zahl, sonst ggf. in den grauen Balken mit Linie zurueck
        if S >= 1.0:
            if S >= MIN_W:
                ax.text(S/2, y[i]+off, f"{S:.1f}", ha="center", va="center", fontsize=5.8)
            else:
                ax.annotate(f"{S:.1f}", xy=(S/2, y[i]+off), xytext=(S+4, y[i]+off),
                            fontsize=5.4, ha="left", va="center",
                            arrowprops=dict(arrowstyle="-", lw=0.5, color="#505050"))
        # B liegt in der Mitte -> bleibt im Balken (in der Praxis nie so schmal)
        if B > 0.8:
            ax.text(S+B/2, y[i]+off, f"{B:.1f}", ha="center", va="center", fontsize=5.8)
        # R liegt rechts im Balken -> bei Bedarf nach rechts auslagern
        if R_ > 0.1:
            if R_ >= MIN_W:
                ax.text(min(S+B+R_/2, 98), y[i]+off, f"{R_:.1f}",
                        ha="center", va="center", fontsize=5.8)
            else:
                ax.annotate(f"{R_:.1f}", xy=(100, y[i]+off), xytext=(103, y[i]+off),
                            fontsize=5.8, ha="left", va="center",
                            arrowprops=dict(arrowstyle="-", lw=0.5, color="#909090"))
ax.set_yticks(y); ax.set_yticklabels([NAME[a] for a in ABS], fontsize=7.5)
ax.set_xlabel("Isolates (%)"); ax.set_xlim(0, 130)
ax.set_xticks([0, 20, 40, 60, 80, 100])
ax.set_title("Report states", loc="left", fontsize=8.5)
ax.legend(handles=[Patch(facecolor=OI["skyblue"], label="non-resistant"),
                   Patch(facecolor=LGREY, label="blank"),
                   Patch(facecolor=OI["vermillion"], label="resistant")],
          fontsize=6.4, handlelength=1.0, ncol=3, loc="upper left",
          bbox_to_anchor=(-0.15, -0.16), columnspacing=1.0)
ax.text(0.0, -0.29, "upper bar: model probability only\n"
        "lower bar: with enzyme rule-in",
        transform=ax.transAxes, ha="left", va="top", fontsize=6.0,
        color="black", linespacing=1.5)
panel_label(ax, "a", dx=-24.0)

# --- b) Beobachtete Resistenz in den berichteten Gruppen --------------------
ax_b = fig.add_subplot(gs[0, 1])
ax = ax_b
for i, ab in enumerate(ABS):
    v = RI[ab]
    if v["Sr"] is not None:
        ax.plot(i-0.20, v["Sr"], "o", color=OI["skyblue"], ms=7)
        ax.text(i-0.20, v["Sr"]+3.5, f"{v['Sr']:.1f}", ha="center", fontsize=6.2)
    if v["Rr_old"] is not None:
        ax.plot(i+0.20, v["Rr_old"], "o", color="white", ms=7,
                mec=OI["vermillion"], mew=1.5)
        ax.text(i+0.20, v["Rr_old"]-3.5, f"{v['Rr_old']:.0f}",
                ha="center", va="top", fontsize=6.2)
    if v["Rr_new"] is not None:
        ax.plot(i+0.20, v["Rr_new"], "o", color=OI["vermillion"], ms=7)
        ax.text(i+0.20, v["Rr_new"]+3.5, f"{v['Rr_new']:.0f}", ha="center", fontsize=6.2)
    if v["Rr_old"] is not None and v["Rr_new"] is not None:
        ax.plot([i+0.20, i+0.20], [v["Rr_old"], v["Rr_new"]],
                color=OI["vermillion"], lw=1.0, alpha=0.6, zorder=0)
ax.axhspan(0, 3, color=OI["skyblue"], alpha=0.16)
ax.axhspan(80, 100, color=OI["vermillion"], alpha=0.12)
ax.set_xticks(range(len(ABS)))
ax.set_xticklabels([NAME[a] for a in ABS], rotation=22, ha="right")
ax.set_xlim(-0.6, len(ABS)-0.4)
ax.set_ylabel("Observed resistance (%)"); ax.set_ylim(0, 100)
ax.set_title("Calibration by report", loc="left", fontsize=8.5)
ax.plot([], [], "o", color=OI["skyblue"], ms=6, label="non-resistant report")
ax.plot([], [], "o", color="white", mec=OI["vermillion"], mew=1.5, ms=6,
        label="resistant report, model only")
ax.plot([], [], "o", color=OI["vermillion"], ms=6,
        label="resistant report, with enzyme rule-in")
hs, ls = ax.get_legend_handles_labels()
ax.legend(handles=hs, fontsize=6.0, handlelength=1.0,
          loc="upper left", bbox_to_anchor=(-0.02, -0.2))
panel_label(ax, "b", dx=-24.0)
pos = ax_b.get_position()
ax_b.set_position([
    pos.x0 - 0.02,
    pos.y0,
    pos.width,
    pos.height
])

# --- c) Die drei Zustaende auf der MIC-Skala ------------------------------
ax = fig.add_subplot(gs[0, 2])
Z = json.load(open(DATA + "drei_zustaende.json"))
SC = {"non-resistant": OI["skyblue"], "blank": GREY, "resistant": OI["vermillion"]}
y = np.arange(len(ABS))[::-1]
for i, ab in enumerate(ABS):
    for k, r in enumerate([r for r in Z[ab] if "median" in r]):
        off = 0.26 - k*0.26
        ax.plot([r["q1"], r["q3"]], [y[i]+off]*2, color=SC[r["state"]], lw=2.2,
                solid_capstyle="butt", alpha=0.55)
        ax.plot(r["median"], y[i]+off, "o", color=SC[r["state"]], ms=5.5,
                mec="white", mew=0.7, zorder=3)
ax.axvline(0, color="k", lw=0.9, ls="--")
ax.set_yticks(y); ax.set_yticklabels([NAME[a] for a in ABS], fontsize=7.5)
ax.set_xlabel("log$_2$ distance from breakpoint")
ax.set_xlim(-6.5, 4.8)
ax.set_xticks([-6, -4, -2, 0, 2, 4])
ax.set_ylim(-0.75, len(ABS)-0.25)
ax.set_title("MIC-scale position", loc="left", fontsize=8.5)
ax.legend(handles=[Patch(facecolor=SC[s_], label=l) for s_, l in
                   [("non-resistant", "reported non-resistant"),
                    ("blank", "left blank"),
                    ("resistant", "reported resistant")]],
          fontsize=6.2, handlelength=1.0, loc="upper left",
          bbox_to_anchor=(-0.02, -0.20))
panel_label(ax, "c", dx=-24.0)

save(fig, OUT + "Figure6_reporting")