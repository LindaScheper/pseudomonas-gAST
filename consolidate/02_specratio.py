"""Antibiotikaspezifitaet eines Markerblocks auf der Odds-Ratio-Skala.

Der bisherige Quotient verglich die positiven praediktiven Werte. Der PPV haengt
aber von der Praevalenz der jeweiligen Substanz ab, und die reicht in dieser
Sammlung von 8.4 % (Tobramycin) bis 19.2 % (Ciprofloxacin). Ein Marker, der die
Resistenzchance fuer alle vier Substanzen um denselben Faktor erhoeht, bekommt auf
der PPV-Skala trotzdem einen Quotienten von 1.6 fuer Ciprofloxacin und 0.6 fuer
Tobramycin. Auf der Odds-Ratio-Skala ist der Quotient in diesem Fall exakt 1.
"""
import pandas as pd, numpy as np, pickle, collections, json
ABS = ["CIP","CAZ","MEM","TOB"]
THR = 1.4

def load():
    C = pickle.load(open("buildcache.pkl","rb"))
    common = set(C["CIP"][0].isolate)
    for ab in ABS: common &= set(C[ab][0].isolate)
    IDX={}; Y={}
    for ab in ABS:
        d,_,_ = C[ab]
        keep=[i for i,iso in enumerate(d.isolate) if iso in common]
        s=d.iloc[keep]; IDX[ab]={iso:k for k,iso in enumerate(s.isolate)}; Y[ab]=s.R.values
    carr=collections.defaultdict(set)
    for ab in ABS:
        d,f2,_=C[ab]
        for i,f in enumerate(f2):
            iso=d.isolate.iloc[i]
            if iso in common:
                for m in f: carr[m].add(iso)
    return C, common, IDX, Y, carr

def contrast(S, IDX, Y):
    """log OR je Substanz, daraus Kontrast eigene Substanz gegen Mittel der drei anderen."""
    L={};V={};N={};P={}
    for ab in ABS:
        y=Y[ab]; m=np.zeros(len(y),bool)
        m[[IDX[ab][x] for x in S if x in IDX[ab]]]=True
        a=y[m].sum()+.5; b=(~y[m]).sum()+.5; c=y[~m].sum()+.5; e=(~y[~m]).sum()+.5
        L[ab]=np.log((a/b)/(c/e)); V[ab]=1/a+1/b+1/c+1/e
        N[ab]=int(m.sum()); P[ab]=float(y[m].mean())
    out={}
    for ab in ABS:
        oth=[o for o in ABS if o!=ab]
        con=L[ab]-np.mean([L[o] for o in oth])
        se=np.sqrt(V[ab]+sum(V[o] for o in oth)/9)   # konservativ: Unabhaengigkeit unterstellt
        out[ab]=dict(n=N[ab], ppv=P[ab], orr=float(np.exp(L[ab])),
                     ratio=float(np.exp(con)), lo=float(np.exp(con-1.96*se)),
                     hi=float(np.exp(con+1.96*se)),
                     ratio_ppv=float(P[ab]/np.mean([P[o] for o in oth])) if np.mean([P[o] for o in oth])>0 else float("inf"))
    return out

def specific(r):
    return r["ratio"] > THR and r["lo"] > 1.0

if __name__ == "__main__":
    C, common, IDX, Y, carr = load()
    PREV={ab: float(Y[ab].mean()) for ab in ABS}
    # Verlaengerungsklasse mitfuehren
    CLS=pickle.load(open("oprd_cls.pkl","rb"))
    carr["ELONG:oprD"]=CLS["alle verlaengert"]&common
    rows=[]
    for mk,S in carr.items():
        if len(S)<20: continue
        r=contrast(S,IDX,Y)
        if min(r[ab]["n"] for ab in ABS)<20: continue
        for ab in ABS:
            rows.append(dict(marker=mk, ab=ab, **r[ab], spec=specific(r[ab])))
    R=pd.DataFrame(rows); R.to_pickle("specratio.pkl")
    # Nullexperiment
    null={}
    for OR in [2,4,8,16,32]:
        p={ab: (PREV[ab]/(1-PREV[ab])*OR)/(1+PREV[ab]/(1-PREV[ab])*OR) for ab in ABS}
        null[OR]={ab: float(p[ab]/np.mean([p[o] for o in ABS if o!=ab])) for ab in ABS}
    J=dict(prev=PREV, n_common=len(common), thr=THR, null=null,
           n_pairs=int(len(R)), n_spec_or=int(R.spec.sum()),
           n_spec_ppv=int((R.ratio_ppv>THR).sum()),
           gained=int((R.spec & (R.ratio_ppv<=THR)).sum()),
           lost=int(((~R.spec) & (R.ratio_ppv>THR)).sum()))
    key=["TRUNC:oprD","LOSS:oprD","ELONG:oprD","gyrA_83","parC_87","TRUNC:PBP4",
         "TRUNC:ampD","fusA1_671","oprJ_25","TRUNC:mexA","TRUNC:mexB","LOSS:mexA","STOP:mexR","STOP:nalD","STOP:nfxB","STOP:ampD","STOP:mucA","FS:galU","STOP:oprD","STOP:mexZ"]
    J["key"]={}
    for mk in key:
        q=R[R.marker==mk]
        if len(q): J["key"][mk]={r.ab: dict(n=int(r.n),ppv=r.ppv,orr=r.orr,ratio=r.ratio,
                                            lo=r.lo,hi=r.hi,ratio_ppv=r.ratio_ppv,spec=bool(r.spec))
                                 for _,r in q.iterrows()}
    json.dump(J, open("specratio.json","w"), indent=1)
    print(f"gemeinsame Isolate {len(common):,}; {len(R)//4} Markerbloecke")
    print(f"spezifisch auf OR-Skala: {J['n_spec_or']}   auf PPV-Skala: {J['n_spec_ppv']}")
    print(f"neu als spezifisch: {J['gained']}   verlieren den Status: {J['lost']}")
    print("\nNullexperiment (Marker mit identischem OR fuer alle vier):")
    for OR,v in null.items():
        print(f"  OR={OR:<3}" + "".join(f"  {ab} {v[ab]:.2f}" for ab in ABS))
    print("\nSchluesselmarker:")
    for mk,v in J["key"].items():
        best=max(v, key=lambda a: v[a]["ratio"])
        r=v[best]
        print(f"  {mk:<14} {best}  n={r['n']:>4}  PPV {r['ppv']*100:5.1f}%  OR {r['orr']:7.2f}  "
              f"Quotient {r['ratio']:6.2f} ({r['lo']:.2f}-{r['hi']:.2f})  [PPV-Skala {r['ratio_ppv']:.2f}]"
              + ("  spezifisch" if r["spec"] else ""))
