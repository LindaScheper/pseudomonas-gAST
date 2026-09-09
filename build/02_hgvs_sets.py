"""Die Mutationstabelle enthaelt HGVS-Notation fuer Frameshift, Stopcodon und Insertion.
Bisher wurden nur einfache Substitutionen und die Proteinlaenge ausgewertet."""
import pandas as pd, numpy as np, re, pickle
M=pd.read_pickle("v3_mut.pkl"); M=M[~M.gene.isin({"mexT","colR"})]
s=M.mut_str.fillna("")
FS  = re.compile(r"([A-Z])(\d+)([A-Z]?)fs")
STOP= re.compile(r"([A-Z])(\d+)\*(?!\d)")     # Y59*  (nicht fs*152)
INS = re.compile(r"(\d+)_?(\d*)ins")
DEL = re.compile(r"([A-Z*])(\d+)-")
for name,rx in [("Frameshift (fs)",FS),("vorzeitiges Stopcodon",STOP),("Insertion (ins)",INS),("Deletion",DEL)]:
    m=s.str.contains(rx)
    print(f"{name:<24}{m.sum():>8,} Isolat-Locus-Paare   {M[m].isolate.nunique():>6,} Isolate   {M[m].gene.nunique():>3} Loci")
print()
ex=s[s.str.contains("ins")].head(3).tolist()
print("Beispiele Insertion:", [e[:120] for e in ex])
# Traegermengen je Locus
def carriers(rx):
    m=s.str.contains(rx)
    return {g:set(x.isolate) for g,x in M[m].groupby("gene")}
CFS=carriers(FS); CST=carriers(STOP); CIN=carriers(INS)
pickle.dump({"fs":CFS,"stop":CST,"ins":CIN}, open("hgvs_sets.pkl","wb"))
print("\nLoci mit den meisten Frameshift-Traegern:")
for g,v in sorted(CFS.items(), key=lambda kv:-len(kv[1]))[:15]:
    print(f"  {g:<10}{len(v):>6}")
print("\nLoci mit den meisten Stopcodon-Traegern:")
for g,v in sorted(CST.items(), key=lambda kv:-len(kv[1]))[:15]:
    print(f"  {g:<10}{len(v):>6}")
