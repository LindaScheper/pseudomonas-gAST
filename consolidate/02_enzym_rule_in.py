"""Der enzymbasierte Rule-in-Pfad, gerechnet auf dem aktuellen Modellauf.

Die Rule-in-Liste ist vorab aus der Substratspezifitaet der Enzymklassen
festgelegt, nicht aus den eigenen praediktiven Werten abgeleitet:

  Meropenem   Metallo-Beta-Laktamasen (IMP, VIM, NDM, FIM, SPM, GIM) und
              Klasse-A-Carbapenemasen (GES-5, GES-13, KPC)
  Ceftazidim  dieselben, dazu die Extended-Spectrum-Enzyme (PER, VEB, PME,
              BEL, GES-1, SHV-2 und die OXA-10-Gruppe mit erweitertem
              Spektrum: OXA-10, -14, -15, -17, -35, -141)
  Tobramycin  16S-rRNA-Methyltransferasen (ArmA, RmtB) und die
              Aminoglykosid-modifizierenden Enzyme mit Tobramycin im
              Substratspektrum: AAC(6')-I, ANT(4')-II, ANT(2'')-I
  Ciprofloxacin  keine

Wichtig: die Wahrscheinlichkeiten werden aus daten/final_all.json gelesen und
nicht neu gerechnet, damit die Zahlen im Manuskript aus genau demselben Lauf
stammen wie alle uebrigen. Erzeugt daten/rulein_v2.json.
"""
import sys, os, json, pickle
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
from paths import DATA

C = pickle.load(open(DATA + "buildcache.pkl", "rb"))
P = json.load(open(DATA + "final_all.json"))

MBL   = ("blaIMP", "blaVIM", "blaNDM", "blaFIM", "blaSPM", "blaGIM")
CARBA = MBL + ("blaGES-5", "blaGES-13", "blaKPC")
ESBL  = ("blaPER", "blaVEB", "blaPME", "blaBEL", "blaGES-1", "blaOXA-10",
         "blaOXA-14", "blaOXA-15", "blaOXA-17", "blaOXA-35", "blaOXA-141", "blaSHV-2")
RULEIN = {"MEM": CARBA, "CAZ": CARBA + ESBL,
          "TOB": ("armA", "rmtB", "aac(6')-I", "ant(4')-II", "ant(2'')-I"),
          "CIP": ()}


def is_rulein(m, ab):
    if not m.startswith("ACQ:"):
        return False
    g = m.replace("ACQ:", "")
    return any(g.startswith(p) for p in RULEIN[ab])


OUT = {}
for ab in ["CIP", "CAZ", "MEM", "TOB"]:
    d, f2, gn = C[ab]
    y = d.R.values
    x = P[ab]["voll"]
    Pr = np.array(x["Pr"])
    assert list(d.isolate) == x["iso"], "Isolatreihenfolge passt nicht"
    ri = np.array([any(is_rulein(m, ab) for m in f) for f in f2])
    lo, hi = Pr < 0.05, Pr > 0.80
    blank = (~lo) & (~hi)
    R2, S2 = hi | ri, lo & ~ri
    B2 = ~R2 & ~S2
    OUT[ab] = dict(
        ri_n=int(ri.sum()),
        ri_ppv=float(y[ri].mean()) if ri.sum() else None,
        ri_share=float(y[ri].sum() / y.sum()) if ri.sum() else 0.0,
        ri_blank=int((ri & blank).sum()),
        R_old=float(hi.mean()), Rr_old=float(y[hi].mean()) if hi.sum() else None,
        R_new=float(R2.mean()), Rr_new=float(y[R2].mean()) if R2.sum() else None,
        S=float(S2.mean()), S_n=int(S2.sum()),
        Sr=float(y[S2].mean()) if S2.sum() else None, B=float(B2.mean()))
    o = OUT[ab]
    print("%-4s Liste n=%3d  davon R %s  Anteil aller R %.1f %%  vom Modell blank %d" % (
        ab, o["ri_n"],
        ("%.1f %%" % (100 * o["ri_ppv"])) if o["ri_ppv"] is not None else "     -",
        100 * o["ri_share"], o["ri_blank"]))
    print("     resistent berichtet %.1f %% (R %.1f %%)  ->  %.1f %% (R %.1f %%)" % (
        100 * o["R_old"], 100 * (o["Rr_old"] or 0), 100 * o["R_new"], 100 * (o["Rr_new"] or 0)))
    print("     nicht-resistent %.1f %% (n=%d, R %.1f %%),  blank %.1f %%" % (
        100 * o["S"], o["S_n"], 100 * (o["Sr"] or 0), 100 * o["B"]))

json.dump(OUT, open(DATA + "rulein_v2.json", "w"), indent=1)
print("\ngeschrieben:", DATA + "rulein_v2.json")
