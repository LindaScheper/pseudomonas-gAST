"""Was kostet der Spezifitaetsfilter an echten Determinanten?

Der Filter wird im Manuskript an einem Gegenbeispiel gezeigt (nalD). Diese
Rechnung beziffert ihn systematisch: von allen Bloecken, die das
Determinantenkriterium erfuellen und an einem der 33 mechanistisch etablierten
Loci des jeweiligen Antibiotikums liegen, wie viele ueberstehen das
Spezifitaetskriterium?

Getrennt ausgewiesen wird, woran ein verworfener Block scheitert:
  - am Abstand   (Quotient <= 1.4, also kaum Praeferenz fuer eine Substanz)
  - an der Schaerfe (Quotient > 1.4, aber untere Konfidenzgrenze <= 1, also
    zu wenige Traeger, um die Praeferenz abzusichern)

Die Locusliste ist dieselbe wie in vorstufen/06_katalogvergleich_33loci.py.

Laeuft direkt aus daten/.  -> stdout und daten/filterkosten.json
"""
import os, sys, json, pickle
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
from paths import DATA

CORE = {"CIP": {"gyrA", "gyrB", "parC", "parE", "mexR", "nalC", "nalD", "nfxB", "mexS", "mexZ"},
        "CAZ": {"ampC", "ampD", "ampDh2", "ampDh3", "ampR", "ampG", "mpl", "nagZ",
                "PBP4", "PBP3", "mexR", "nalC", "nalD", "oprM"},
        "MEM": {"oprD", "ampC", "ampD", "ampR", "mexR", "nalC", "nalD", "mexS"},
        "TOB": {"mexZ", "armZ", "mexX", "mexY", "fusA1", "parR", "parS", "pmrA",
                "pmrB", "amgS", "phoP", "phoQ"}}

S = pickle.load(open(DATA + "specratio.pkl", "rb"))


def locus(m):
    m = m.split(" [")[0]
    for p in ("TRUNC:", "LOSS:", "FS:", "STOP:", "ELONG:"):
        if m.startswith(p):
            return m[len(p):]
    if m.startswith("ACQ:"):
        return "ACQ"
    return m.rsplit("_", 1)[0]


print("Bloecke an mechanistisch etablierten Loci, die das Determinantenkriterium")
print("erfuellen (mindestens zehn Traeger, mindestens die Haelfte resistent).\n")
print("%-14s %8s %9s %11s %12s" % ("", "Bloecke", "behalten", "am Abstand", "an der Schaerfe"))
OUT = {}
gk = gn = ga = gs = 0
for ab in ["CIP", "CAZ", "MEM", "TOB"]:
    sub = S[(S.ab == ab) & (S.n >= 10) & (S.ppv >= 0.5)]
    sub = sub[[locus(m) in CORE[ab] for m in sub.marker]]
    n = len(sub)
    k = int(sub.spec.sum())
    ab_ = int(((~sub.spec) & (sub.ratio <= 1.4)).sum())
    sc = int(((~sub.spec) & (sub.ratio > 1.4)).sum())
    OUT[ab] = dict(n=n, behalten=k, am_abstand=ab_, an_der_schaerfe=sc)
    gn += n; gk += k; ga += ab_; gs += sc
    print("%-14s %8d %6d %2.0f %% %6d %2.0f %% %7d %2.0f %%" % (
        ab, n, k, 100*k/n, ab_, 100*ab_/n, sc, 100*sc/n))
OUT["gesamt"] = dict(n=gn, behalten=gk, am_abstand=ga, an_der_schaerfe=gs)
print("%-14s %8d %6d %2.0f %% %6d %2.0f %% %7d %2.0f %%" % (
    "gesamt", gn, gk, 100*gk/gn, ga, 100*ga/gn, gs, 100*gs/gn))

json.dump(OUT, open(DATA + "filterkosten.json", "w"), indent=1)
print("\ngeschrieben:", DATA + "filterkosten.json")
print("\nAn der Schaerfe scheitern heisst: die Praeferenz fuer eine Substanz ist da,")
print("aber der Block hat zu wenige Traeger, um sie abzusichern. Diese Bloecke")
print("gingen bei einer groesseren Sammlung wieder in das Panel ein.")
