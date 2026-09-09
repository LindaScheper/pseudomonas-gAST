"""Figure 3: Vorhersagbarkeit folgt dem Mechanismus.
Eine Tafel je Antibiotikum, identische Achsen: praediktiver Wert einer
Mechanismusklasse gegen den Anteil der Resistenz, den sie erfasst.
Symbolflaeche proportional zur Traegerzahl.
Die Aussage liegt in der oberen rechten Ecke: eine Klasse, die praezise ist
und zugleich breit abdeckt. Bei Ciprofloxacin liegt dort das Wirkziel; bei
den drei uebrigen Substanzen ist die Ecke leer.
Loest die frueher auf Ceftazidim beschraenkte Fassung ab.
Eingabe : mech_classes.json (vorstufen/15_mechanismusklassen.py), bss aus figdata.pkl
Ausgabe : Figure3_mechanism_classes.png / .pdf
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import *
from paths import DATA, OUT
import json, pickle
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
ABS  = ["CIP", "CAZ", "MEM", "TOB"]
NAME = {"CIP": "Ciprofloxacin", "CAZ": "Ceftazidime",
        "MEM": "Meropenem", "TOB": "Tobramycin"}
CLS  = {"Acquired enzymes":  "#000000",
        "Drug target":       OI["vermillion"],
        "AmpC pathway":      OI["green"],
        "Efflux regulation": OI["blue"],
        "Efflux structural": OI["skyblue"],
        "Porin / LPS":       OI["orange"],
        "Cell wall, other":  GREY}
SHORT = {"Acquired enzymes": "enzymes", "Drug target": "target",
         "AmpC pathway": "ampC", "Efflux regulation": "efflux reg.",
         "Efflux structural": "efflux genes", "Porin / LPS": "porin",
         "Cell wall, other": "other"}
M   = json.load(open(DATA + "mech_classes.json"))
BSS = pickle.load(open(DATA + "figdata.pkl", "rb"))["bss"]

# Manueller Versatz je Panel (Antibiotikum) und Mechanismusklasse, in Punkten
# (dx, dy) relativ zum Kreismittelpunkt. Startwerte grob nach Augenmass -
# bitte pro Panel einzeln nachjustieren, bis nichts mehr ueberlappt.
# Fehlt ein Eintrag, wird DEFAULT_OFFSET verwendet.
DEFAULT_OFFSET = (0, 10)
LABEL_OFFSET = {
    ("CIP", "Drug target"):        (0, 14),
    ("CIP", "Efflux structural"):  (0, 10),
    ("CIP", "Efflux regulation"):  (2, -14),

    ("CAZ", "Acquired enzymes"):   (0, 10),
    ("CAZ", "Drug target"):        (0, 10),
    ("CAZ", "AmpC pathway"):       (17, 0),
    ("CAZ", "Efflux regulation"):  (14, 7),
    ("CAZ", "Porin / LPS"):        (4, 7),
    ("CAZ", "Efflux structural"):  (16, -9),
    ("CAZ", "Cell wall, other"):   (0, -10),

    ("MEM", "Drug target"):        (0, 10),
    ("MEM", "Acquired enzymes"):   (19, 2),
    ("MEM", "AmpC pathway"):       (0, -11),
    ("MEM", "Efflux regulation"):  (0, 12),
    ("MEM", "Efflux structural"):  (15, -8),
    ("MEM", "Porin / LPS"):        (0, -12),

    ("TOB", "Porin / LPS"):        (0, 9),
    ("TOB", "Acquired enzymes"):   (20, 4),
    ("TOB", "Drug target"):        (0, 8),
    ("TOB", "Efflux regulation"):  (0, 10),
    ("TOB", "Efflux structural"):  (0, -10),
}

fig = plt.figure(figsize=(15*CM, 13.5*CM))
gs  = GridSpec(2, 2, figure=fig, hspace=0.46, wspace=0.30)
XLIM, YLIM = (0, 82), (33, 105)
for k, ab in enumerate(ABS):
    ax = fig.add_subplot(gs[k // 2, k % 2])
    ax.set_xlim(*XLIM); ax.set_ylim(*YLIM)
    # Zielbereich: praezise und abdeckend zugleich
    ax.add_patch(plt.Rectangle((50, 80), XLIM[1]-50, 100-80,
                           facecolor=OI["yellow"], alpha=0.16, zorder=0))
    for j, r in enumerate(sorted(M[ab], key=lambda r: r["cover"])):
        x, y = r["cover"]*100, r["ppv"]*100
        rad = 14 + r["carr"]*0.22
        ax.scatter(x, y, s=rad, color=CLS[r["lab"]],
                   alpha=0.82, edgecolor="white", linewidth=0.7, zorder=3)
        dx, dy = LABEL_OFFSET.get((ab, r["lab"]), DEFAULT_OFFSET)
        ax.annotate(SHORT[r["lab"]], (x, y), textcoords="offset points",
                    xytext=(dx, dy), ha="center", va="center",
                    fontsize=6.0, color=CLS[r["lab"]], zorder=4)
    ax.set_yticks([40,60,80,100])
    ax.axhline(80, color=LGREY, lw=0.7, zorder=0)
    ax.vlines(50, YLIM[0], 100, color=LGREY, lw=0.7, zorder=0)
    ax.set_title(NAME[ab], loc="left", fontsize=8.2)
    if k % 2 == 0:
        ax.set_ylabel("Resistant Carriers (%)")
    if k // 2 == 1:
        ax.set_xlabel("Share of Resistance covered (%)")
    panel_label(ax, "abcd"[k], dx=-24)

# Groessenlegende: zeigt, welche Traegerzahl welcher Symbolgroesse entspricht
size_vals = [50, 200, 800]
size_handles = [Line2D([], [], marker="o", ls="", color=GREY, markeredgecolor="white",
                       markersize=np.sqrt(14 + v*0.22), label=str(v))
                for v in size_vals]

# Klassen-Legende: eine Farbe je Mechanismusklasse
h = [Line2D([], [], marker="o", ls="", color=CLS[c], ms=5, label=c) for c in CLS]

# Groessenlegende (oben)
leg1 = fig.legend(handles=size_handles, loc="center left", bbox_to_anchor=(1.02, 0.70),
                   fontsize=6.6, handlelength=1.6, labelspacing=1.3,
                   frameon=False, title="Carriers (n)", title_fontsize=6.6,
                   alignment="left")

# Klassen-Legende (darunter)
leg2 = fig.legend(handles=h, loc="center left", bbox_to_anchor=(1.02, 0.45),
                   fontsize=6.6, handlelength=1.0, labelspacing=1.0, frameon=False,
                   title="Mechanism Class", title_fontsize=6.6,
                   alignment="left")
# Platz rechts fuer die Legenden schaffen
fig.subplots_adjust(right=1.0)
save(fig, OUT + "Figure3_mechanism_classes")