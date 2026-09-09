"""Isolate mit gestoerter MexXY-Pumpe: Trägerzahl, Resistenz, Kontrollrechnungen.

Liefert die Zahlen des MexXY-Absatzes im Ergebnisteil. Die Trägerdefinition ist
dieselbe wie in vorstufen/pumps.py — Trunkierung, Genverlust, Frameshift oder
Stopcodon in mexX oder mexY — hier aber aus buildcache.pkl statt aus den
Rohdaten, damit die Rechnung ohne die Rohtabellen nachvollziehbar bleibt. Die
Trägerzahl reproduziert pumps.json exakt (362).

Laeuft direkt aus daten/.  -> stdout
"""
import os, sys, pickle, json
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
from paths import DATA

C  = pickle.load(open(DATA + "buildcache.pkl", "rb"))
FA = json.load(open(DATA + "final_all.json"))
d, f2, _ = C["TOB"]
y = dict(zip(d.isolate, d.R.values.astype(bool)))

VERLUST = ("TRUNC:", "LOSS:", "FS:", "STOP:")
AME = {"ACQ:rmtB4", "ACQ:aac(6')-Ib", "ACQ:aac(3)-IIa", "ACQ:aac(3)-IIe",
       "ACQ:ant(2'')-Ia", "ACQ:ant(4')-IIb", "ACQ:armA"}


def gestoert(gene):
    s = set()
    for iso, fs in zip(d.isolate, f2):
        for m in fs:
            for p in VERLUST:
                if m.startswith(p) and m[len(p):].split(" [")[0] in gene:
                    s.add(iso); break
            else:
                continue
            break
    return s


def quote(s):
    s = list(s)
    return len(s), 100 * sum(y[i] for i in s) / len(s)


def odds(a, b):
    na, ra = len(a), sum(y[i] for i in a)
    nb, rb = len(b), sum(y[i] for i in b)
    return (ra / (na - ra)) / (rb / (nb - rb))


xy   = gestoert({"mexX", "mexY"})
rest = set(d.isolate) - xy
ame  = {iso for iso, fs in zip(d.isolate, f2) if fs & AME}

print("MexXY gestoert : n = %d, %.1f %% resistent" % quote(xy))
print("uebrige        : n = %d, %.1f %% resistent" % quote(rest))
print("Odds Ratio     : %.2f  (ungeschichtet; geschichtet 0.33, siehe pumps.json)" % odds(xy, rest))

print("\nDavon Traeger einer erworbenen aminoglykosidmodifizierenden Enzymklasse: %d" % len(xy & ame))
a, b = xy - ame, rest - ame
print("ohne diese     : gestoert n = %d, %.1f %% | uebrige n = %d, %.1f %% | Odds Ratio %.2f"
      % (quote(a) + quote(b) + (odds(a, b),)))

pat = {i: str(p) for i, p in zip(d.isolate, d.patient)}
st  = {i: str(s) for i, s in zip(d.isolate, d.st)}
P = {pat[i] for i in xy} - {"nan", "None", ""}
S = {st[i] for i in xy} - {"-", "NA", "nan"}
print("\nTraeger verteilt auf %d Patienten und %d Sequenztypen" % (len(P), len(S)))

for rule, name in [("dedup_first", "erstes Isolat"), ("dedup_worst", "resistentestes")]:
    keep = set(FA["TOB"][rule]["iso"])
    c = [i for i in xy if i in keep]
    print("nach Deduplikation (%s): n = %d, %.1f %% resistent"
          % (name, len(c), 100 * sum(y[i] for i in c) / len(c)))
