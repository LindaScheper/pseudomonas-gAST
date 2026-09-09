"""Abbildung 5: MIC-abhaengige Detektion, Sockel und Decke.

Eingabe : figdata.pkl
Ausgabe : Figure5_mic_detection.png / .pdf
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

# =========================== Abbildung 4 ===========================
fig = plt.figure(figsize=(20*CM, 8.5*CM))
gs = GridSpec(1, 3, figure=fig, wspace=0.52, width_ratios=[1.25, 1, 1])

ax = fig.add_subplot(gs[0, 0])
for ab in ["CIP", "CAZ", "MEM", "TOB"]:
    m = D["detect"][ab]["mic"]
    if not m: continue
    bp = D["detect"][ab]["bp"]
    xs = [np.log2(v/bp) for v, n, d, r in m]
    ys = [d*100 for v, n, d, r in m]
    ax.plot(xs, ys, "o-", color=COL[ab], ms=4, lw=1.4, label=NAME[ab])
ax.axvline(0, color="k", lw=0.8, ls="--")
ax.set_xlabel("log$_2$ distance from breakpoint")
ax.set_ylabel("Isolates carrying a determinant (%)")
ax.set_ylim(0, 110); ax.legend(loc="upper left", fontsize=7, handlelength=1.4)
panel_label(ax, "a", dx=-24.0)

ax = fig.add_subplot(gs[0, 1])
m = D["detect"]["CIP"]["canon"]
bp = D["detect"]["CIP"]["bp"]
xs_mic = [np.log2(v) for v, n, d in m]
ax.plot(xs_mic, [d*100 for v, n, d in m],
        "o-", color=OI["vermillion"], ms=4.5, lw=1.5)
ax.axvline(np.log2(bp), color="k", lw=0.8, ls="--")
ax.set_xticks(xs_mic)
ax.set_xticklabels([f"{v:g}" for v, n, d in m], rotation=35, ha="right", fontsize=7)

# NEU: Tick-Labels um 6 pt nach rechts verschieben
import matplotlib.transforms as mtransforms
offset = mtransforms.ScaledTranslation(4/72, 0, fig.dpi_scale_trans)
for label in ax.get_xticklabels():
    label.set_transform(label.get_transform() + offset)

ax.set_xlabel("MIC (mg/L)")
ax.set_ylabel(r"Isolates carrying $\it{gyrA}$/$\it{parC}$" "\n" r"target mutation (%)")
ax.set_ylim(0, 110)
panel_label(ax, "b", dx=-24.0)

ax = fig.add_subplot(gs[0, 2])

labs, fp, fn, abs_present = [], [], [], []

for ab in ["CIP", "CAZ", "MEM", "TOB"]:
    m = D["detect"][ab]["mic"]
    if not m:
        continue

    labs.append(NAME[ab])
    abs_present.append(ab)
    fp.append(D["detect"][ab]["floor"] * 100)
    fn.append((1 - m[-1][2]) * 100)

x = np.arange(len(labs))
width = 0.36


# False positive:
# gestrichelt/schraffiert in Antibiotikafarbe,
# aber mit durchgezogener Außenlinie
for i, ab in enumerate(abs_present):
    ax.bar(
        x[i] - 0.19,
        fp[i],
        width,
        facecolor="white",
        edgecolor=COL[ab],
        linewidth=1.2,
        hatch="//"
    )


# False negative:
# vollständig mit Antibiotikafarbe gefüllt
for i, ab in enumerate(abs_present):
    ax.bar(
        x[i] + 0.19,
        fn[i],
        width,
        facecolor=COL[ab],
        edgecolor=COL[ab],
        linewidth=1.0
    )


# Werte über den Balken
for i, (a, b) in enumerate(zip(fp, fn)):
    ax.text(
        i - 0.19,
        a + 0.8,
        f"{a:.0f}",
        ha="center",
        fontsize=6.8
    )

    ax.text(
        i + 0.19,
        b + 0.8,
        f"{b:.0f}",
        ha="center",
        fontsize=6.8
    )


ax.set_xticks(x)
ax.set_xticklabels(labs, rotation=25, ha="right")

# NEU: Tick-Labels um 4 pt nach rechts verschieben
offset_c = mtransforms.ScaledTranslation(4/72, 0, fig.dpi_scale_trans)
for label in ax.get_xticklabels():
    label.set_transform(label.get_transform() + offset_c)

ax.set_ylabel("Isolates misclassified by the marker test (%)")
ax.set_ylim(0, 35)


# Legende
from matplotlib.patches import Patch

ax.legend(
    handles=[
        Patch(
            facecolor="white",
            edgecolor="black",
            linewidth=1.0,
            hatch="//////",
            label="false positive"
        ),
        Patch(
            facecolor="black",
            edgecolor="black",
            label="false negative"
        )
    ],
    loc="upper left",
    fontsize=6.2,
    handlelength=1.4,
    labelspacing=0.4,
    frameon=False,
    bbox_to_anchor=(-0.03, 1.02)
)

panel_label(ax, "c", dx=-24.0)
save(fig, OUT + "Figure5_mic_detection")
