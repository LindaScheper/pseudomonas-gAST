"""Effluxpumpen: Welche Stoerung senkt, welche erhoeht die Resistenz je Substanz?

Ausgeschlossen sind die 63 Isolate, bei denen die Annotation weitgehend ausgefallen
ist; sie wuerden fuer jede Pumpe einen falschen Verlust beitragen.
"""
import pickle, numpy as np, pandas as pd, collections, json
from scipy.stats import fisher_exact
C=pickle.load(open("buildcache.pkl","rb")); H=pickle.load(open("hgvs_sets.pkl","rb"))
M=pd.read_pickle("v3_mut.pkl")
genes=sorted(M.gene.unique()); present={g:set(s.isolate) for g,s in M.groupby("gene")}
nmiss={i:sum(1 for g in genes if i not in present[g]) for i in sorted(set(M.isolate))}
BAD={i for i,n in nmiss.items() if n>=20}
MODE=M.dropna(subset=["plen"]).groupby("gene").plen.agg(lambda s:s.mode().iloc[0])
Mp=M.dropna(subset=["plen"]).copy(); Mp["rel"]=Mp.plen/Mp.gene.map(MODE)
TRUNC={g:set(s.isolate) for g,s in Mp[Mp.rel<0.95].groupby("gene")}
ABSENT={g:(set(M.isolate)-present[g]) for g in genes}
def disrupted(gs):
    s=set()
    for g in gs:
        s|=TRUNC.get(g,set()) | ABSENT.get(g,set()) | H["fs"].get(g,set()) | H["stop"].get(g,set())
    return s - BAD
PUMPS={
 "MexAB-OprM": ["mexA","mexB"],
 "MexCD-OprJ": ["mexC","mexD","oprJ"],
 "MexEF-OprN": ["mexE","mexF","oprN"],
 "MexXY":      ["mexX","mexY"],
}
NAME={"CIP":"Ciprofloxacin","CAZ":"Ceftazidime","MEM":"Meropenem","TOB":"Tobramycin"}
OUT={}
print("Odds Ratio fuer Resistenz bei gestoerter Pumpe, Mantel-Haenszel ueber Sequenztypen")
print("unter 1 = Stoerung macht empfindlicher, ueber 1 = Stoerung geht mit mehr Resistenz einher\n")
print(f"{'Pumpe':<14}{'n':>6}" + "".join(f"{NAME[a]:>16}" for a in ["CIP","CAZ","MEM","TOB"]))
print("-"*84)
for pump,gs in PUMPS.items():
    S=disrupted(gs); row=[]; rec={}
    for ab in ["CIP","CAZ","MEM","TOB"]:
        d,f2,gn=C[ab]; y=d.R.values
        keep=[i for i,s in enumerate(d.isolate) if s not in BAD]
        d2=d.iloc[keep]; y2=y[keep]
        car=np.array([s in S for s in d2.isolate])
        if car.sum()<20: row.append("  zu wenige"); continue
        st=d2.st.astype(str).values
        num=den=0.0
        for s_ in pd.unique(st):
            j=st==s_
            A=int((car&j&y2).sum()); B=int((car&j&~y2).sum())
            Cc=int((~car&j&y2).sum()); Dd=int((~car&j&~y2).sum())
            n=A+B+Cc+Dd
            if n and (A+B) and (Cc+Dd): num+=A*Dd/n; den+=B*Cc/n
        mh=num/den if den>0 else np.nan
        a=int((car&y2).sum()); b=int(car.sum())-a
        c=int((~car&y2).sum()); dd=int((~car).sum())-c
        _,p=fisher_exact([[a,b],[c,dd]])
        star="*" if p<0.01 else " "
        row.append(f"{mh:>8.2f}{star} ({a/max(a+b,1)*100:4.1f}%)")
        rec[ab]=dict(n=int(car.sum()), ppv=a/max(a+b,1), mh=float(mh), p=float(p))
    print(f"{pump:<14}{len(S):>6}" + "".join(f"{x:>16}" for x in row))
    OUT[pump]=rec
json.dump(OUT,open("pumps.json","w"),indent=1)
print("\n* p < 0.01 (Fisher, ungeschichtet).  In Klammern der Anteil resistenter Traeger.")
print(f"\nAusgeschlossen: {len(BAD)} Isolate mit >= 20 fehlenden Loci")
