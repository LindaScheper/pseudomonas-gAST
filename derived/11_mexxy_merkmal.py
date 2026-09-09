"""Tobramycin mit einem vorzeichenbehafteten Merkmal fuer die Stoerung von MexXY.

Gegenstueck zu 10_pumpenmerkmal.py: dort MexAB-OprM fuer Meropenem, hier MexXY
fuer Tobramycin. Der Vergleichslauf ohne das Merkmal steht in final_all.json
unter TOB/voll.

Im Arbeitsverzeichnis mit den Rohdaten-Pickles aufrufen.  -> tob_xyprot.json
"""
import sys, json, pickle
sys.path.insert(0, ".")
import dedup_run as D

C = pickle.load(open("buildcache.pkl", "rb"))
D.build = lambda ab: (C[ab][0], C[ab][1], {}, C[ab][2])
XY = ("TRUNC:mexX", "LOSS:mexX", "FS:mexX", "STOP:mexX",
      "TRUNC:mexY", "LOSS:mexY", "FS:mexY", "STOP:mexY")
D.PROTECT = lambda ff: 1.0 if any(m in ff for m in XY) else 0.0

x = D.run("TOB", True, 3, None)
json.dump({"TOB": x}, open("tob_xyprot.json", "w"), indent=1)
print("TOB mit MexXY-Merkmal: BSS %.3f  rule-out %.1f %% (R %.1f %%)" % (
    x["bss"], 100 * x["lo"], 100 * (x["lo_r"] or 0)))
