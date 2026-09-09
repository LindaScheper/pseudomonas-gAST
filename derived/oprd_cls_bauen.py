"""Erzeugt oprd_cls.pkl fuer 02_specratio.py.

02_specratio.py erwartet unter dem Schluessel "alle verlaengert" die Menge
der Isolat-IDs, bei denen oprD verlaengert ist (Proteinlaenge > 105 % des
Populationsmodus). Das ist dieselbe Definition, die 07_elongation.py fuer
seine Verlaengerungsklasse verwendet (dort: rel > 1.05, gruppiert nach Gen).

Aufruf im Arbeitsverzeichnis, nach 01_tabelle_einlesen.py / 03_zusammenfuehren.py:
    python3 oprd_cls_bauen.py
"""
import pandas as pd, pickle

M = pd.read_pickle("v3_mut.pkl")
M = M[~M.gene.isin({"mexT", "colR"})]
MODE = M.dropna(subset=["plen"]).groupby("gene").plen.agg(lambda s: s.mode().iloc[0])
M = M.dropna(subset=["plen"]).copy()
M["rel"] = M.plen / M.gene.map(MODE)

long_oprd = set(M[(M.gene == "oprD") & (M.rel > 1.05)].isolate)

CLS = {"alle verlaengert": long_oprd}
pickle.dump(CLS, open("oprd_cls.pkl", "wb"))
print(f"oprd_cls.pkl geschrieben: {len(long_oprd):,} Isolate mit verlaengertem oprD")
