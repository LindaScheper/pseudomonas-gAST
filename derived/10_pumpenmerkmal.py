"""Meropenem ohne das vorzeichenbehaftete Pumpenmerkmal.

Der Lauf mit Merkmal steht in final_all.json unter MEM/voll. Die gepaarte
Differenz der beiden Brier Skill Scores wird in
pruefungen/01_konfidenzintervalle.py gebildet; hier wird nur der Vergleichslauf
erzeugt.

Im Arbeitsverzeichnis mit den Rohdaten-Pickles aufrufen.  -> noprot.json
"""
import sys, json, pickle
sys.path.insert(0, ".")
import dedup_run as D

C = pickle.load(open("buildcache.pkl", "rb"))
D.build = lambda ab: (C[ab][0], C[ab][1], {}, C[ab][2])
D.PROTECT = None

x = D.run("MEM", True, 3, None)
json.dump({"MEM": x}, open("noprot.json", "w"), indent=1)
print("MEM ohne Pumpenmerkmal: BSS %.3f" % x["bss"])
