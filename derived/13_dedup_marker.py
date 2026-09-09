"""Praediktiver Wert einzelner Determinanten vor und nach Deduplikation.

Liefert die Zahlen fuer Supplementary Figure 4c und fuer die Saetze dazu im
Text. Frueher standen diese Werte als feste Liste in 01_figdata.py; sie waren
fuer zwei Marker (trunkiertes PBP4, fusA1 Position 154) nicht mehr aktuell.

Die deduplizierten Isolatmengen werden nicht neu gebildet, sondern aus
final_all.json uebernommen — das sind genau die Mengen, auf denen auch die
Brier Skill Scores im Supplement beruhen.

Laeuft direkt aus daten/.  -> daten/dedup_markers.json
"""
import os, sys, pickle, json
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
from paths import DATA

C  = pickle.load(open(DATA + "buildcache.pkl", "rb"))
FA = json.load(open(DATA + "final_all.json"))

# Marker, Substanz, Praefix des Blocknamens in der Merkmalsmenge
ROWS = [("parC S87",      "CIP", "parC_87"),
        ("gyrA T83",      "CIP", "gyrA_83"),
        ("rmtB4",         "TOB", "ACQ:rmtB4"),
        ("aac(6')-Ib",    "TOB", "ACQ:aac(6')-Ib"),
        ("ampD trunc.",   "CAZ", "TRUNC:ampD"),
        ("PBP4 trunc.",   "CAZ", "TRUNC:PBP4"),
        ("oprD trunc.",   "MEM", "TRUNC:oprD"),
        ("nalD trunc.",   "MEM", "TRUNC:nalD"),
        ("mexZ pos.126",  "TOB", "mexZ_126"),
        ("fusA1 pos.154", "TOB", "fusA1_154")]


def traeger(ab, pref):
    """Isolate, die den Block tragen. Kollabierte Bloecke heissen 'name [+n]'."""
    d, f2, _ = C[ab]
    return {iso for iso, s in zip(d.isolate, f2)
            if any(m == pref or m.startswith(pref + " [") for m in s)}


OUT = []
print("%-15s %-4s | %-16s | %-16s | %-16s" % (
    "Marker", "AB", "volle Sammlung", "erstes Isolat", "resistentestes"))
print("-" * 78)
for lab, ab, pref in ROWS:
    d, _, _ = C[ab]
    y = dict(zip(d.isolate, d.R.values.astype(bool)))
    car = traeger(ab, pref)
    v = {}
    for rule in ["voll", "dedup_first", "dedup_worst"]:
        keep = set(FA[ab][rule]["iso"])
        c = [i for i in car if i in keep]
        v[rule] = (len(c), 100 * sum(y[i] for i in c) / len(c) if c else float("nan"))
    print("%-15s %-4s | %5d  %5.1f %%   | %5d  %5.1f %%   | %5d  %5.1f %%" % (
        lab, ab, v["voll"][0], v["voll"][1], v["dedup_first"][0], v["dedup_first"][1],
        v["dedup_worst"][0], v["dedup_worst"][1]))
    OUT.append([lab, round(v["voll"][1], 1), round(v["dedup_first"][1], 1),
                v["dedup_first"][0], round(v["dedup_worst"][1], 1), v["dedup_worst"][0]])

json.dump(OUT, open(DATA + "dedup_markers.json", "w"), indent=1)
print("\ngeschrieben:", DATA + "dedup_markers.json")
print("Die Abbildung zeigt die erste Regel (chronologisch erstes Isolat).")
