"""Fuehrt die Teildateien aus 01_tabelle_einlesen.py zusammen.

01_tabelle_einlesen.py wird mehrfach aufgerufen (Zeitlimit pro Lauf) und
schreibt dabei pro Lauf eigene Teildateien:
    v3_part_<n>.pkl   Mutationstabelle (DataFrame: gene, isolate, mut_str, plen, reflen)
    v3_ph_<n>.pkl     Phaenotyp-Dict je Isolat (v_AB, n_AB, bp_AB, t_AB fuer AB in
                       MEM, CAZ, CIP, TOB, GEN)

Dieses Skript fasst alle vorhandenen Teile zusammen zu:
    v3_mut.pkl     DataFrame, wie von 03_core.py / 04_buildcache.py erwartet
    v3_pheno.pkl   DataFrame mit Spalte "isolate" + v_/n_/bp_/t_<AB>-Spalten

Aufruf erst NACHDEM 01_tabelle_einlesen.py "alle Blaetter gelesen" gemeldet hat.
"""
import glob, pickle, pandas as pd, json, os

STATE = "v3_state.json"
if os.path.exists(STATE):
    st = json.load(open(STATE))
    print(f"Stand laut {STATE}: {len(st.get('done', []))} Blaetter erledigt")

part_files = sorted(glob.glob("v3_part_*.pkl"))
ph_files   = sorted(glob.glob("v3_ph_*.pkl"))
if not part_files or not ph_files:
    raise SystemExit("Keine v3_part_*.pkl / v3_ph_*.pkl gefunden - erst 01_tabelle_einlesen.py laufen lassen.")

print(f"{len(part_files)} Mutations-Teildateien, {len(ph_files)} Phaenotyp-Teildateien gefunden")

# --- Mutationstabelle: DataFrames aneinanderhaengen ---
M = pd.concat([pd.read_pickle(f) for f in part_files], ignore_index=True)
M = M.drop_duplicates(subset=["gene", "isolate", "mut_str"])
M.to_pickle("v3_mut.pkl")
print(f"v3_mut.pkl geschrieben: {len(M):,} Zeilen, {M.isolate.nunique():,} Isolate, {M.gene.nunique()} Gene")

# --- Phaenotyp: Dicts zusammenfuehren (je Isolat gewinnt der erste Treffer, wie im Originalskript) ---
pheno = {}
for f in ph_files:
    d = pickle.load(open(f, "rb"))
    for iso, vals in d.items():
        if iso not in pheno:
            pheno[iso] = vals
P = pd.DataFrame.from_dict(pheno, orient="index")
P.index.name = "isolate"
P = P.reset_index()
P.to_pickle("v3_pheno.pkl")
print(f"v3_pheno.pkl geschrieben: {len(P):,} Isolate")
