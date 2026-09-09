"""Die drei Befundzustaende sind geordnet — der Blank ist kein Nichtwissen.

Fuer jede Substanz werden die drei Zustaende des Befunds, die sich aus der
Modellwahrscheinlichkeit ergeben (unter 5 %, dazwischen, ueber 80 %), gegen
den Phaenotyp gestellt:

  - beobachtete Resistenz je Zustand, ueber die ganze Sammlung
  - Lage auf der MIC-Skala, nur Zentrum A, wo abgestufte MIC-Werte vorliegen:
    Median und Quartile des vorzeichenbehafteten log2-Abstands vom Grenzwert,
    und der Anteil der Isolate im Fenster von einer Verduennungsstufe um den
    Grenzwert

Die Blank-Gruppe liegt in beidem zwischen den beiden berichteten Gruppen und
ist fuer grenzwertnahe Isolate angereichert. Der dritte Zustand sagt also
nicht "keine Aussage moeglich", sondern "dieser Stamm liegt im Uebergang".

Laeuft direkt aus daten/.  -> daten/drei_zustaende.json
"""
import os, sys, json, pickle
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
from paths import DATA

C = pickle.load(open(DATA + "buildcache.pkl", "rb"))
F = json.load(open(DATA + "final_all.json"))
ABS = ["CIP", "MEM", "CAZ", "TOB"]
NAME = {"CIP": "Ciprofloxacin", "MEM": "Meropenem",
        "CAZ": "Ceftazidime", "TOB": "Tobramycin"}
STATES = [("non-resistant", lambda p: p < 0.05),
          ("blank",         lambda p: (p >= 0.05) & (p <= 0.80)),
          ("resistant",     lambda p: p > 0.80)]

OUT = {}
for ab in ABS:
    d = C[ab][0]
    A = F[ab]["voll"]
    iso = A["iso"]; p = np.array(A["Pr"])
    ymap = dict(zip(d.isolate, d.R.values.astype(float)))
    nmap = dict(zip(d.isolate, d["n_%s" % ab]))
    cmap = dict(zip(d.isolate, d.centre))
    y = np.array([ymap[i] for i in iso])
    # Zentrum A: abgestufte MIC-Werte
    selA = np.array([cmap[i] == "A" and nmap[i] == nmap[i] for i in iso])
    norm = np.array([nmap[i] if nmap[i] == nmap[i] else np.nan for i in iso])

    rows = []
    print("\n%s   (%d Isolate, Praevalenz %.1f %%)" % (NAME[ab], len(y), 100 * y.mean()))
    print("   %-14s %6s %8s | %s" % ("Zustand", "n", "R", "Zentrum A: Median  Q1..Q3   im Fenster +-1"))
    for lab, fn in STATES:
        m = fn(p)
        r = dict(state=lab, n=int(m.sum()), share=float(m.mean()),
                 resistant=float(y[m].mean()) if m.sum() else None)
        mA = m & selA
        if mA.sum() >= 5:
            v = norm[mA]
            r.update(nA=int(mA.sum()), median=float(np.median(v)),
                     q1=float(np.percentile(v, 25)), q3=float(np.percentile(v, 75)),
                     near=float(np.mean(np.abs(v) <= 1)))
            print("   %-14s %6d %7.1f%% | %+6.1f  %+5.1f..%+5.1f  %5.1f %%  (n=%d)"
                  % (lab, r["n"], 100 * r["resistant"], r["median"], r["q1"], r["q3"],
                     100 * r["near"], r["nA"]))
        else:
            r.update(nA=int(mA.sum()))
            print("   %-14s %6d %7.1f%% | zu wenige Isolate mit MIC" % (lab, r["n"], 100 * r["resistant"]))
        rows.append(r)
    OUT[ab] = rows

json.dump(OUT, open(DATA + "drei_zustaende.json", "w"), indent=1)
print("\ngeschrieben:", DATA + "drei_zustaende.json")
print("\nDie Zustaende beruhen auf der Modellwahrscheinlichkeit allein; der")
print("Enzym-Rule-in verschiebt zusaetzlich Isolate aus dem Blank in den")
print("Resistenzbefund und wuerde die Ordnung nur verstaerken.")
