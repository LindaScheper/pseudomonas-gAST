"""Wieviel bringt der Antibiotika-Spezifitaetsfilter?

Wiederholt die vollstaendige Kreuzvalidierung einmal ohne den Filter, damit sie
gegen den Lauf mit Filter (final_all.json, Eintrag "voll") gehalten werden kann.

Wie die uebrigen Skripte in diesem Ordner im Arbeitsverzeichnis mit den
Rohdaten-Pickles aufrufen. Ein Lauf je Aufruf, weil ein Lauf rund 90 Sekunden
braucht; so oft aufrufen, bis "nichts offen" erscheint.

  -> spec_off.json
"""
import sys, os, json, pickle
sys.path.insert(0, ".")
import dedup_run as D

C = pickle.load(open("buildcache.pkl", "rb"))
D.build = lambda ab: (C[ab][0], C[ab][1], {}, C[ab][2])

PORIN = ("TRUNC:oprD", "LOSS:oprD")
PUMP  = ("TRUNC:mexA", "LOSS:mexA", "TRUNC:mexB", "LOSS:mexB")
D.PROTECT = lambda ff: 1.0 if (any(m in ff for m in PUMP)
                               and not any(m in ff for m in PORIN)) else 0.0

F = "spec_off.json"
R = json.load(open(F)) if os.path.exists(F) else {}
todo = [ab for ab in ["CIP", "CAZ", "MEM", "TOB"] if ab not in R]
if not todo:
    print("nichts offen"); raise SystemExit

ab = todo[0]
x = D.run(ab, False, 3, None)          # use_spec=False
R[ab] = x
json.dump(R, open(F, "w"), indent=1)
print("%s ohne Filter: BSS %.3f  k50 %.1f/%.1f   [%d offen]" % (
    ab, x["bss"], x["fr"][50]["sens"] * 100, x["fr"][50]["spec"] * 100, len(todo) - 1))
