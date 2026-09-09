"""Abbildung: Effluxpumpen — Substratmatrix, Porindominanz, Modellgewinn."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import *
from paths import DATA, OUT
import json, numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Rectangle

J=json.load(open(DATA+"pumpfig.json"))

# Tafel c kommt aus dem aktuellen Lauf, nicht aus dem aelteren model-Block in
# pumpfig.json: ohne Merkmal aus noprot.json, mit Merkmal aus final_all.json.
# Das sind dieselben Zahlen, die pruefungen/01_konfidenzintervalle.py verwendet
# und die im Manuskript stehen (0.213 -> 0.223, 38.4 % -> 41.9 %).
def _mem_modell():
    import pickle
    C=pickle.load(open(DATA+"buildcache.pkl","rb"))
    d=C["MEM"][0]; ymap=dict(zip(d.isolate, d.R.values.astype(float)))
    A=json.load(open(DATA+"final_all.json"))["MEM"]["voll"]
    B=json.load(open(DATA+"noprot.json"))["MEM"]
    assert A["iso"]==B["iso"], "unterschiedliche Isolatreihenfolge"
    y=np.array([ymap[i] for i in A["iso"]], float)
    out={}
    for tag,v in [("Basismodell",B),("+ mexA/mexB gestoert, Porin intakt",A)]:
        p=np.array(v["Pr"]); s=p<0.05
        out[tag]=dict(bss=1-np.mean((p-y)**2)/np.mean((y.mean()-y)**2),
                      lo=float(s.mean()), lo_r=float(y[s].mean()))
    return out
J["model"]=_mem_modell()
ABS=["CIP","CAZ","MEM","TOB"]
NAME={"CIP":"Ciprofloxacin","CAZ":"Ceftazidime","MEM":"Meropenem","TOB":"Tobramycin"}
COL={"CIP":OI["vermillion"],"CAZ":OI["green"],"MEM":OI["blue"],"TOB":OI["purple"]}
PUMPS=["MexAB-OprM","MexCD-OprJ","MexEF-OprN","MexXY"]

# manuell vorgegebene Zuordnung Pumpe -> zugehörige(s) Antibiotikum/Antibiotika
PUMP_ABX = {
    "MexAB-OprM": ["MEM", "CAZ", "CIP"],
    "MexCD-OprJ": ["CIP"],
    "MexEF-OprN": ["CIP"],
    "MexXY":      ["TOB"],
}

fig=plt.figure(figsize=(18*CM, 14*CM))
gs=GridSpec(2,2,figure=fig,hspace=0.55,wspace=0.62,height_ratios=[1,0.85],width_ratios=[1.05,1])

# ---------- a) Matrix
ax=fig.add_subplot(gs[0,0])
Mx=np.array([[J["matrix"][p][a]["mh"] for a in ABS] for p in PUMPS])
im=ax.imshow(np.log2(Mx), cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-2, vcenter=0, vmax=2.4), aspect="auto")
# signifikant nach Bonferroni-Korrektur ueber alle 16 Kombinationen
Pv=np.array([[J["matrix"][p][a]["p"] for a in ABS] for p in PUMPS])
for i in range(len(PUMPS)):
    for j in range(len(ABS)):
        v=Mx[i,j]
        sig="*" if Pv[i,j] < 0.05/16 else ""
        ax.text(j,i,f"{v:.2f}{sig}",ha="center",va="center",fontsize=7.6,
                color="white" if abs(np.log2(v))>1.2 else "#202020",
                fontweight="normal", zorder=6)
ax.set_xticks(range(4)); ax.set_xticklabels([NAME[a] for a in ABS],rotation=22,ha="right",fontsize=7.5)
ax.set_yticks(range(4)); ax.set_yticklabels(PUMPS,fontsize=7.2)
ax.set_title("Odds Ratio for Resistance when a Pump is disrupted",loc="left",fontsize=8.5)

# weiße Umrandung der Zellen, die zur jeweiligen Pumpe gehören
for i, p in enumerate(PUMPS):
    for a in PUMP_ABX[p]:
        j = ABS.index(a)
        rect = Rectangle((j-0.5, i-0.5), 1, 1,
                          fill=False, edgecolor="black",
                          linewidth=1.8, zorder=5)
        rect.set_clip_on(False)
        ax.add_patch(rect)

cb=fig.colorbar(im,ax=ax,orientation="horizontal",fraction=0.055,pad=0.3,
                ticks=np.log2([0.25,0.5,1,2,4]))
cb.ax.set_xticklabels(["0.25","0.5","1","2","4"],fontsize=6.8)
cb.set_label("Odds Ratio for resistance, stratified by sequence type", fontsize=6.4)
panel_label(ax,"a",dx=-24)

# ---------- b) Befundgewinn
ax=fig.add_subplot(gs[1,0])
base=J["model"]["Basismodell"]; new=J["model"]["+ mexA/mexB gestoert, Porin intakt"]
x=np.arange(2); w=0.34
ax.bar(x-w/2,[base["lo"]*100,base["bss"]*100],w,color=GREY,label="markers of resistance only")
ax.bar(x+w/2,[new["lo"]*100,new["bss"]*100],w,color=OI["blue"],label="with the protective marker")
for i,(a,b) in enumerate([(base["lo"]*100,new["lo"]*100),(base["bss"]*100,new["bss"]*100)]):
    ax.text(i-w/2,a+0.8,f"{a:.1f}",ha="center",fontsize=6.8)
    ax.text(i+w/2,b+0.8,f"{b:.1f}",ha="center",fontsize=6.8)
ax.set_xticks(x); ax.set_xticklabels(["reportable as\nnon-resistant (%)","Brier skill\nscore (×100)"],fontsize=7.2)
ax.set_ylim(0,52)
ax.set_ylim(0,56)
ax.set_title("Marker's effect on meropenem report",loc="left",fontsize=8.5)
panel_label(ax,"b",dx=-24)

# ---------- c) MIC-Modell
ax=fig.add_subplot(gs[1,1])
mb=J["mic"]["Basismodell"]; mn=J["mic"]["mit Marker"]
ax.bar([0-w/2],[mb["r2"]],w,color=GREY,label="markers of resistance only")
ax.bar([0+w/2],[mn["r2"]],w,color=OI["blue"],label="with the protective marker")
ax.text(0-w/2,mb["r2"]+0.008,f"{mb['r2']:.2f}",ha="center",fontsize=7.2)
ax.text(0+w/2,mn["r2"]+0.008,f"{mn['r2']:.2f}",ha="center",fontsize=7.2)
ax.set_xlim(-0.5,0.5); ax.set_ylim(0,0.40)
ax.set_xticks([])
ax.set_ylabel("Explained Variance (R²)")
ax.set_title("Marker's effect on predicted MIC",loc="left",fontsize=8.5)
panel_label(ax,"c",dx=-24)
from matplotlib.patches import Patch
fig.legend(handles=[Patch(facecolor=GREY,label="resistance marker only"),
                    Patch(facecolor=OI["blue"],label="+ pump/porin feature")],
           loc="lower center", bbox_to_anchor=(0.5,0.0), ncol=2, fontsize=7,
           handlelength=1.2, columnspacing=2.0, frameon=False)
save(fig,OUT+"Figure4_pumps")