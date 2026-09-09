"""Woher kommt die 1.4 — und wie viel haengt an ihr?

Die Schwelle stammt aus der Zeit, als der Spezifitaetsquotient auf der
Skala der praediktiven Werte gebildet wurde. Dort gibt ein Marker OHNE jede
Spezifitaet — also mit identischem Odds Ratio fuer alle vier Substanzen —
allein wegen der unterschiedlichen Praevalenzen einen Quotienten bis 1.59
zurueck. 1.4 war der Abstand, den ein voellig unspezifischer Marker auf
dieser Skala noch erreichen kann.

Auf der Odds-Ratio-Skala, die das Manuskript jetzt verwendet, faellt diese
Begruendung weg: dort gibt derselbe Marker exakt 1.00 zurueck, unabhaengig
von der Praevalenz. Die 1.4 ist damit eine reine Konvention.

Dieses Skript zeigt beides: das Nullexperiment, aus dem die Zahl stammt, und
wie wenig am konkreten Wert haengt.

Laeuft direkt aus daten/.  -> stdout
"""
import os, sys, json, pickle
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
from paths import DATA

J = json.load(open(DATA + "specratio.json"))
S = pickle.load(open(DATA + "specratio.pkl", "rb"))
ABS = ["CIP", "CAZ", "MEM", "TOB"]

print("Nullexperiment: Quotient eines Markers mit identischem Odds Ratio fuer alle vier,")
print("gebildet auf der Skala der PRAEDIKTIVEN WERTE. Praevalenzen %s\n"
      % {k: round(v * 100, 1) for k, v in J["prev"].items()})
print("%-6s" % "OR" + "".join("%8s" % a for a in ABS))
for OR, v in J["null"].items():
    print("%-6s" % OR + "".join("%8.2f" % v[a] for a in ABS))
print("\nEin voellig unspezifischer Marker erreicht so bis zu %.2f (Ciprofloxacin, OR 2)."
      % max(v["CIP"] for v in J["null"].values()))
print("Daher die Schwelle 1.4. Auf der Odds-Ratio-Skala waere derselbe Wert exakt 1.00.\n")

print("Wie viel haengt am konkreten Wert?")
print("Kriterium: Quotient > Schwelle UND untere 95-%-Grenze > 1\n")
print("%9s %12s %10s" % ("Schwelle", "spezifisch", "Anteil"))
for thr in [1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 2.0, 3.0]:
    k = int(((S.ratio > thr) & (S.lo > 1)).sum())
    print("%9.1f %12d %9.1f %%" % (thr, k, 100 * k / len(S)))

sub = S[S.lo > 1]
print("\nBloecke mit unterer Grenze > 1, nach Quotient:")
for a, b in [(1.0, 1.2), (1.2, 1.4), (1.4, 1.6), (1.6, 2.0), (2.0, 1e9)]:
    n = int(((sub.ratio >= a) & (sub.ratio < b)).sum())
    print("   %4.1f bis %-5s %5d" % (a, ("%.1f" % b) if b < 1e9 else "hoeher", n))
print("\nZwischen 1.0 und 1.2 liegt kein einziger Block: ein Quotient dicht an 1 kann")
print("kein Konfidenzintervall haben, das die 1 ausschliesst. Die Konfidenzgrenze")
print("erledigt die Arbeit, der Abstand entscheidet nur ueber %d von %d Bloecken."
      % (int(((S.ratio > 1.0) & (S.lo > 1)).sum()) - int(((S.ratio > 1.4) & (S.lo > 1)).sum()),
         int(((S.ratio > 1.0) & (S.lo > 1)).sum())))
