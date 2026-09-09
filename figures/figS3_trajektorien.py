"""Supplementary Figure 3: Sensitivitaet gegen 1 minus Spezifitaet, Panelgroesse
10 bis 120 Bloecke, kreuzvalidiert.
Zeigt, dass die Reihenfolge der vier Substanzen an jedem Betriebspunkt dieselbe
bleibt. Der im Manuskript berichtete Punkt liegt bei 50 Bloecken.
War frueher Tafel a von Abbildung 4.
Eingabe : figdata.pkl
Ausgabe : SupplementaryFigure3_trajectories.png / .pdf
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import *
from paths import DATA, OUT
import pickle
import numpy as np
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
COL  = {"CIP": OI["vermillion"], "CAZ": OI["green"], "MEM": OI["blue"], "TOB": OI["purple"]}
NAME = {"CIP": "Ciprofloxacin", "CAZ": "Ceftazidime", "MEM": "Meropenem", "TOB": "Tobramycin"}
ABS  = ["CIP", "CAZ", "MEM", "TOB"]
KS   = [10, 20, 30, 50, 80, 120]
D = pickle.load(open(DATA + "figdata.pkl", "rb"))
fig, ax = plt.subplots(figsize=(11*CM, 8.5*CM))

def shade(color, s):
    """s=0 -> weiss, s=1 -> Grundfarbe, s=2 -> schwarz; alpha bleibt bei 1 (deckend)"""
    rgb = np.array(mcolors.to_rgb(color))
    if s <= 1:
        white = np.array([1.0, 1.0, 1.0])
        return tuple(white + s * (rgb - white))
    else:
        return tuple(rgb * (2 - s))

def darken(color, frac=0.6):
    """frac=1 -> unveraendert, frac=0 -> schwarz; alpha bleibt bei 1 (deckend)"""
    rgb = np.array(mcolors.to_rgb(color))
    return tuple(rgb * frac)

n_ks = len(KS)
t = np.linspace(0, 1, n_ks) ** 1.6
shades = 0.2 + 1.2 * t   # 0.2 (blass) ... 1.4 (dunkler als die Grundfarbe)
for ab in ABS:
    fr = D["frontier"][ab]["mit"]["fr"]
    xs = [(1 - fr[str(k)]["spec"]) * 100 for k in KS]
    ys = [fr[str(k)]["sens"] * 100 for k in KS]
    ax.plot(xs, ys, "-", color=COL[ab], lw=1.6, label=NAME[ab], zorder=1)
    for x, y, s in zip(xs, ys, shades):
        ax.plot(x, y, "o", color=shade(COL[ab], s),
                markeredgecolor=darken(COL[ab], 0.45), markeredgewidth=0.8,
                ms=4.5, zorder=3)

ax.plot([0, 100], [0, 100], color=LGREY, lw=0.8, ls="--", zorder=0)
ax.set_xlabel("1 − specificity (%)"); ax.set_ylabel("Sensitivity (%)")
ax.set_xlim(0, 72); ax.set_ylim(0, 100)

# Legende 1: Antibiotika (Farbe = Substanz)
leg1 = ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0),
                  fontsize=7.5, handlelength=1.4, frameon=False)
ax.add_artist(leg1)  # damit die zweite Legende diese nicht ueberschreibt

# Legende 2: Blockgroesse (Farbintensitaet = Panelgroesse), neutral in Grau
intensity_handles = [Line2D([0], [0], marker="o", linestyle="none",
                             markerfacecolor=shade(GREY, s),
                             markeredgecolor=darken(GREY, 0.45),
                             markeredgewidth=0.8, markersize=4.5,
                             label=str(k))
                      for k, s in zip(KS, shades)]
ax.legend(handles=intensity_handles, title="Panel size",
          loc="upper left", bbox_to_anchor=(1.02, 0.55),
          fontsize=6.3, title_fontsize=6.8,
          handlelength=1.0, labelspacing=0.35, frameon=False)

fig.subplots_adjust(left=0.14, right=0.72, bottom=0.15)
save(fig, OUT + "SupplementaryFigure3_trajectories")