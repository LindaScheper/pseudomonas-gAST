"""Meropenem: klinischer Breakpoint gegen den epidemiologischen Cut-off.

norm = log2(MIC / 8), also
    klinisch resistent (MIC > 8 mg/L)          <=>  norm >  0
    nicht Wildtyp      (MIC > 2 mg/L, ECOFF)   <=>  norm > -2

Nur Zentrum A, weil dort graduierte MIC-Werte vorliegen; fuer Hemmhof-
durchmesser gibt es keinen entsprechenden Schnitt. Alles uebrige — Katalog,
Kollabierung, Spezifitaetsfilter, Faltung, Merkmale — bleibt unveraendert,
getauscht wird ausschliesslich der vorhergesagte Endpunkt.

Im Arbeitsverzeichnis mit den Rohdaten-Pickles aufrufen. Ein Lauf je Aufruf;
zweimal aufrufen, bis "nichts offen" erscheint.   -> ecoff_mem.json
"""
import sys, json, pickle, os
sys.path.insert(0, ".")
import dedup_run as D

C = pickle.load(open("buildcache.pkl", "rb"))
d0, f20, gn0 = C["MEM"]
mask = (d0.type == "mic_V").values
d0 = d0[mask].reset_index(drop=True)
f20 = [f for f, m in zip(f20, mask) if m]

PORIN = ("TRUNC:oprD", "LOSS:oprD")
PUMP  = ("TRUNC:mexA", "LOSS:mexA", "TRUNC:mexB", "LOSS:mexB")
D.PROTECT = lambda ff: 1.0 if (any(m in ff for m in PUMP)
                               and not any(m in ff for m in PORIN)) else 0.0

F = "ecoff_mem.json"
R = json.load(open(F)) if os.path.exists(F) else {}
todo = [(t, c) for t, c in [("klinisch", 0.0), ("ecoff", -2.0)] if t not in R]
if not todo:
    print("nichts offen"); raise SystemExit

tag, cut = todo[0]
d = d0.copy()
d["R"] = d.norm > cut
D.build = lambda ab: (d, f20, {}, gn0)
x = D.run("MEM", True, 3, None)
R[tag] = x
json.dump(R, open(F, "w"), indent=1)
print("%-9s n=%d  Praevalenz %.1f %%  BSS %.3f  k50 %.1f/%.1f  "
      "nicht-resistent %.1f %% (R %.1f %%)   [%d offen]" % (
          tag, x["n"], 100 * x["prev"], x["bss"],
          x["fr"][50]["sens"] * 100, x["fr"][50]["spec"] * 100,
          100 * x["lo"], 100 * (x["lo_r"] or 0), len(todo) - 1))
