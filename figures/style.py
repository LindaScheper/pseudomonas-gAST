"""Shared style for the P. aeruginosa gAST mock figures."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Okabe-Ito colourblind-safe palette
OI = dict(black="#000000", orange="#E69F00", skyblue="#56B4E9", green="#009E73",
          yellow="#F0E442", blue="#0072B2", vermillion="#D55E00", purple="#CC79A7")
GREY = "#9A9A9A"; LGREY = "#DCDCDC"

DRUGS = ["Meropenem", "Ciprofloxacin", "Ceftazidime", "Tobramycin"]
DRUG_C = {"Meropenem": OI["blue"], "Ciprofloxacin": OI["vermillion"],
          "Ceftazidime": OI["green"], "Tobramycin": OI["purple"]}

# Marker-load colours (0 / 1 / >=2)
LOAD_C = [OI["skyblue"], OI["orange"], OI["vermillion"]]

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8,
    "axes.labelsize": 8.5, "axes.titlesize": 9.5, "axes.titleweight": "bold",
    "axes.titlepad": 8.0,   # neu: einheitlicher Abstand Titel <-> Plot (in Punkten)
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "figure.dpi": 110, "savefig.dpi": 300, "savefig.bbox": "tight",
    "legend.frameon": False,
})

CM = 1 / 2.54
from matplotlib.transforms import ScaledTranslation
def panel_label(ax, letter, dx=-24):
    fig = ax.figure
    fig.canvas.draw()  # noetig, damit die Titel-Bbox verfuegbar ist
    title = ax.title
    bbox = title.get_window_extent(renderer=fig.canvas.get_renderer())
    bbox_axes = bbox.transformed(fig.transFigure.inverted())
    # Titel-Oberkante in Figure-Koordinaten -> zurueck in Punkte relativ zur Achse
    ax_bbox = ax.get_position()
    dy_pt = (bbox_axes.y1 - ax_bbox.y1) * fig.get_size_inches()[1] * 72
    trans = ax.transAxes + ScaledTranslation(dx/72, dy_pt/72, ax.figure.dpi_scale_trans)
    ax.text(0, 1, letter, transform=trans, fontsize=11,
            fontweight="bold", va="baseline", ha="left")

def mock_banner(fig, text="MOCK DATA - placeholder values for layout review"):
    fig.text(0.995, 0.002, text, ha="right", va="bottom", fontsize=6.5,
             color=OI["vermillion"], style="italic", alpha=0.9)

def save(fig, name):
    fig.savefig(f"{name}.png"); fig.savefig(f"{name}.pdf")
    plt.close(fig); print(f"  wrote {name}.png / .pdf")
