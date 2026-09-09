"""Neue Datenbasis: 95 Gene, oprD beide Zentren, AMRFinder komplett, norm_value."""
import sys; sys.path.insert(0,".")
import pandas as pd, numpy as np, collections, re, pickle
from collapse import collapse_features
M=pd.read_pickle("v3_mut.pkl"); P=pd.read_pickle("v3_pheno.pkl"); A=pd.read_pickle("v3_amr.pkl")
S=pd.read_pickle("sample_table.pkl"); PA=S[S.species=="P. aeruginosa"].drop_duplicates("sample")
PA["sample"]=PA["sample"].astype(str)
META=PA.set_index("sample")[["patient","mlst","date","hospital","Num_contigs","coverage"]].to_dict("index")

# In der Extraktion vom 24. August 2026 sind alle 95 Loci auswertbar. colR und
# mexT waren zuvor durch einen Alignment-Kollaps unbrauchbar (colR: 1.109 Aufrufe
# im Median fuer ein 228-Reste-Protein, Positionen bis 1046; mexT: 40 Aufrufe fuer
# 305 Reste). Nach der Neuextraktion liegen beide bei 0 Aufrufen im Median und
# keiner Position jenseits der Proteinlaenge.
# --- Qualitaetsausschluss (Pruefung durch L. Falgenhauer, 29. August 2026)
# 63 Isolate, die die Genom-QC nicht bestanden haben oder deren Laborspezies nicht
# zur genomischen Vorhersage passt. Sie sind ueber die Probentabelle in die Analyse
# gelangt, weil dort die Laborangabe steht und nicht das QC-Urteil.
#   22  Assemblierungs-QC nicht bestanden (Contigzahl oder Genomgroesse)
#   22  nie durch die WGS-QC-Pipeline gelaufen
#   17  Laborspezies E. coli oder K. pneumoniae, genomisch P. aeruginosa
#    2  Laborspezies P. aeruginosa, genomisch P. paraeruginosa
import json as _json, os as _os
QC_EXCLUDE = set(_json.load(open("exclude63.json"))) if _os.path.exists("exclude63.json") else set()

BROKEN=set()
M=M[~M.gene.isin(BROKEN)].copy() if BROKEN else M.copy()
M["muts"]=M.mut_str.apply(lambda s:[x.strip() for x in str(s).split(",") if x.strip()])
POS=re.compile(r'(\d+)')
# --- 1) Umpolarisierung gegen das Mehrheitsallel je (Gen,Position) ---
cov=M.groupby("gene").isolate.nunique().to_dict()
alt=collections.defaultdict(collections.Counter)
for g,lst in zip(M.gene,M.muts):
    for m in lst:
        d=POS.search(m)
        if d: alt[(g,int(d.group(1)))][m]+=1
major={}
for (g,p),c in alt.items():
    n_alt=sum(c.values()); ref_n=cov[g]-n_alt
    b,bn=c.most_common(1)[0]
    major[(g,p)]=("REF",ref_n) if ref_n>=bn else (b,bn)
flip={k for k,v in major.items() if v[0]!="REF"}
FLIPG=collections.defaultdict(list)
for (g,p) in flip: FLIPG[g].append(p)
DER=collections.defaultdict(set)
for g,iso,lst in zip(M.gene,M.isolate,M.muts):
    seen={}
    for m in lst:
        d=POS.search(m)
        if d: seen[int(d.group(1))]=m
    for p,st in seen.items():
        if major[(g,p)][0]!=st: DER[iso].add(f"{g}_{p}")
    for p in FLIPG.get(g,()):
        if p not in seen: DER[iso].add(f"{g}_{p}")
# --- 2) Trunkierung gegen den Populationsmodus ---
MODE=M.dropna(subset=["plen"]).groupby("gene").plen.agg(lambda s:s.mode().iloc[0])
M["rel"]=M.plen/M.gene.map(MODE)
for g,iso in zip(M[M.rel<0.95].gene,M[M.rel<0.95].isolate): DER[iso].add(f"TRUNC:{g}")
PRESENT={g:set(s) for g,s in M.groupby("gene").isolate.apply(set).items()}
pickle.dump((dict(DER),PRESENT,len(flip),len(major)),open("v2_der.pkl","wb"))
# --- 3) Erworbene Gene aus AMRFinder, substanzspezifisch ---
TOBOK=re.compile(r"^(aac\(6'\)|aac\(3\)-(II|III|IV|VI)|aac\(2'\)|ant\(2''\)|ant\(4'\)|aph\(2''\)|rmt|armA|npmA)")
def acq_for(ab):
    d=collections.defaultdict(set)
    a=A[A.typ=="AMR"]
    if ab in ("CAZ","MEM"):
        # AMRFinder fuehrt bei P. aeruginosa ALLE Beta-Laktamasen als scope="core",
        # weil PDC (chromosomale AmpC) und die OXA-50-Familie intrinsisch sind.
        # Ein Filter auf scope=="plus" liefert deshalb null Treffer. Stattdessen
        # werden intrinsische Familien und Punktmutationen explizit ausgeschlossen.
        b=a[a.cls.astype(str).str.contains("BETA-LACTAM",na=False)].copy()
        b["gene"]=b.gene.astype(str)
        b=b[b.gene.str.startswith("bla")]
        # Klassifikation ueber den Element name, nicht ueber die Genbezeichnung:
        # "OXA-50 family" und die PDC-/class-C-Familie sind bei P. aeruginosa intrinsisch.
        en=b["Element name"].astype(str)
        intrinsic=en.str.contains("OXA-50",na=False)|en.str.contains("PDC",na=False)|\
                  en.str.contains("class C beta-lactamase",na=False)
        sel=b[~intrinsic]
    elif ab=="CIP":
        sel=a[a.cls.astype(str).str.contains("QUINOLONE",na=False)&(a.scope=="plus")]
    else:
        sel=a[a.cls.astype(str).str.contains("AMINOGLYCOSIDE",na=False)&
              a.gene.astype(str).str.match(TOBOK)]
    for iso,g in zip(sel.isolate,sel.gene): d[iso].add("ACQ:"+str(g))
    return d
PANEL={
 "CIP":['gyrA','gyrB','parC','parE','mexR','nalC','nalD','nfxB','mexS','mexC','mexD','oprJ',
        'mexE','mexF','oprN','mexT','mexZ','armZ','oprM','mexA','mexB'],
 "CAZ":['ampC','ampD','ampDh2','ampDh3','ampR','ampG','ampO','ampP','mpl','nagZ','nagA','anmK',
        'dacC','creB','creC','creD','PBP1A','PBP1B','PBP2','PBP3','PBP3A','PBP4','PBP7',
        'mexR','nalC','nalD','oprM','mexA','mexB','galU','oprD','ldcA','LdtPae1','LdtPae2','LdtPae3',
        'slt','sltB1','sltB2','sltB3','mltA','mltD','mltF','mltG','amiA','amiB','amiC','nlpD'],
 "MEM":['oprD','ampC','ampD','ampDh2','ampDh3','ampR','ampG','mpl','nagZ','PBP1A','PBP2','PBP3',
        'mexR','nalC','nalD','oprM','mexA','mexB','mexE','mexF','oprN','mexS','mexT','galU','creB','creC','PBP4'],
 "TOB":['mexZ','armZ','mexX','mexY','oprM','fusA1','parR','parS','pmrA','pmrB','galU',
        'amgR','amgS','amgK','phoP','phoQ','colS','colR'],
}
# --- Frameshift- und Stopcodon-Aufrufe (HGVS-Notation in mut_str)
import re as _re, pickle as _pk, os as _os
if _os.path.exists("hgvs_sets.pkl"):
    _H=_pk.load(open("hgvs_sets.pkl","rb")); _FS=_H["fs"]; _STOP=_H["stop"]
else:
    _FS={}; _STOP={}

def build(ab, min_cov_frac=0.0):
    genes=[g for g in PANEL[ab] if g in PRESENT]
    gs=set(genes); acq=acq_for(ab)
    # Plausibilitaetsfilter: Hemmhofdurchmesser < 6 mm sind unmoeglich (Blaettchen = 6 mm).
    # Betroffen sind 22 Datensaetze, in denen offenbar ein MIC-Wert in der Hemmhofspalte steht;
    # sie wuerden als extrem resistent gewertet, sind aber meist empfindlich.
    P.loc[(P[f"t_{ab}"]=="zone_V")&(P[f"v_{ab}"]<6),[f"n_{ab}",f"v_{ab}"]]=np.nan
    d=P[P[f"v_{ab}"].notna()&P[f"n_{ab}"].notna()].copy()
    d["isolate"]=d.isolate.astype(str)
    d=d[d.isolate.isin(set(M.isolate))]
    d=d[~d.isolate.isin(QC_EXCLUDE)]
    d["norm"]=d[f"n_{ab}"]; d["type"]=d[f"t_{ab}"]
    d["centre"]=np.where(d["type"]=="mic_V","A","B")
    d["R"]=d.norm>0                                  # norm>0 == oberhalb des Breakpoints
    d["st"]=[str(META.get(i,{}).get("mlst","NA")) for i in d.isolate]
    d["patient"]=[str(META.get(i,{}).get("patient","")) for i in d.isolate]
    feat=[]
    for i in d.isolate:
        f={m for m in DER.get(i,set())
           if (m.split(":")[-1] if m.startswith("TRUNC:") else m.rsplit("_",1)[0]) in gs}
        for g in genes:
            if i not in PRESENT.get(g,set()): f.add(f"LOSS:{g}")
        # Frameshift und vorzeitiges Stopcodon aus der HGVS-Notation der Mutationstabelle.
        # Beides inaktiviert ein Protein unabhaengig davon, ob es dabei kuerzer oder
        # laenger wird, und ist an der Proteinlaenge allein nicht zu erkennen.
        for g in genes:
            if i in _FS.get(g,()):   f.add(f"FS:{g}")
            if i in _STOP.get(g,()): f.add(f"STOP:{g}")
        f|=acq.get(i,set()); feat.append(f)
    f2,blocks=collapse_features(feat,0.90,5)
    return d.reset_index(drop=True), f2, blocks, genes
