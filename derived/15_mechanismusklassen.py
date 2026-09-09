"""Mechanismusklassen fuer alle vier Antibiotika.

Bisher gab es diese Auswertung nur fuer Ceftazidim (D["caz_classes"] in
01_figdata.py). Fuer die Thesenabbildung wird sie hier auf alle vier
Substanzen ausgedehnt, mit einer gemeinsamen Klasseneinteilung, damit die
vier Tafeln vergleichbar sind.

Je Klasse werden ausgewiesen:
  n      Zahl der Determinantenbloecke in der Klasse
  carr   Zahl der Isolate, die mindestens einen davon tragen
  ppv    Anteil dieser Traeger, die resistent sind
  cover  Anteil aller resistenten Isolate, die die Klasse erfasst

Determinantenkriterium wie im Manuskript: mindestens zehn Traeger und
mindestens die Haelfte resistent; erworbene Enzyme zaehlen unabhaengig von
der Traegerzahl.

Laeuft direkt aus daten/.  -> daten/mech_classes.json
"""
import os, sys, json, pickle, collections
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
from paths import DATA

C = pickle.load(open(DATA + "buildcache.pkl", "rb"))
ABS = ["CIP", "CAZ", "MEM", "TOB"]

# gemeinsame Klasseneinteilung, mechanisch begruendet
TARGET  = {"gyrA", "gyrB", "parC", "parE", "fusA1",
           "PBP1A", "PBP1B", "PBP2", "PBP3", "PBP3A", "PBP4", "PBP7", "dacC"}
AMPC    = {"ampC", "ampD", "ampDh2", "ampDh3", "ampR", "ampG", "ampO", "ampP",
           "mpl", "nagZ", "nagA", "anmK", "creB", "creC", "creD"}
EFF_REG = {"mexR", "nalC", "nalD", "nfxB", "mexS", "mexT", "mexZ", "armZ",
           "parR", "parS", "pmrA", "pmrB", "amgR", "amgS", "amgK",
           "phoP", "phoQ", "colR", "colS"}
EFF_STR = {"mexA", "mexB", "mexC", "mexD", "oprJ", "mexE", "mexF", "oprN",
           "mexX", "mexY", "oprM"}
PORIN   = {"oprD", "galU"}

ORDER = [("Acquired enzymes", None),
         ("Drug target",      TARGET),
         ("AmpC pathway",     AMPC),
         ("Efflux regulation", EFF_REG),
         ("Efflux structural", EFF_STR),
         ("Porin / LPS",      PORIN),
         ("Cell wall, other", None)]


def gene_of(m):
    if m.startswith("ACQ:"):
        return "ACQ"
    m = m.split(" [")[0]
    for p in ("TRUNC:", "LOSS:", "FS:", "STOP:", "ELONG:"):
        if m.startswith(p):
            return m[len(p):]
    return m.rsplit("_", 1)[0]


def determinanten(f2, y, min_n=10, min_ppv=0.5):
    car = collections.defaultdict(list)
    for i, f in enumerate(f2):
        for m in f:
            car[m].append(i)
    DET = {m for m, c in car.items() if len(c) >= min_n and y[c].mean() >= min_ppv}
    DET |= {m for m in car if m.startswith("ACQ:")}
    return DET, car


OUT = {}
for ab in ABS:
    d, f2, gn = C[ab]
    y = d.R.values
    R = int(y.sum())
    DET, car = determinanten(f2, y)
    rest = set(gn) - TARGET - AMPC - EFF_REG - EFF_STR - PORIN
    rows = []
    print("\n%s  (%d resistente Isolate, %d Determinantenbloecke)" % (ab, R, len(DET)))
    for lab, S in ORDER:
        if lab == "Acquired enzymes":
            ms = [m for m in DET if gene_of(m) == "ACQ"]
        elif lab == "Cell wall, other":
            ms = [m for m in DET if gene_of(m) in rest]
        else:
            ms = [m for m in DET if gene_of(m) in S]
        cc = set()
        for m in ms:
            cc |= set(car[m])
        cc = sorted(cc)
        if not cc:
            continue
        rows.append(dict(lab=lab, n=len(ms), carr=len(cc),
                         ppv=float(y[cc].mean()), cover=float(y[cc].sum() / R)))
        print("   %-19s Bloecke %3d  Traeger %5d  PPV %5.1f %%  Abdeckung %5.1f %%"
              % (lab, len(ms), len(cc), 100 * y[cc].mean(), 100 * y[cc].sum() / R))
    OUT[ab] = rows

json.dump(OUT, open(DATA + "mech_classes.json", "w"), indent=1)
print("\ngeschrieben:", DATA + "mech_classes.json")
print("\nDie Klassen ueberlappen: ein Isolat kann Bloecke aus mehreren Klassen")
print("tragen. Abdeckung und praediktiver Wert sind deshalb je Klasse zu lesen,")
print("nicht zu addieren.")
