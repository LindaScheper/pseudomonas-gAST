import sys; sys.path.insert(0,".")
exec(open("03_core.py").read())
import numpy as np, collections, pickle, json
from sklearn.linear_model import LogisticRegression
from dedup import keep_set
ABS=["CIP","CAZ","MEM","TOB"]
ALL=sorted(set(P.isolate.astype(str))); IX={s:i for i,s in enumerate(ALL)}
HAS=np.zeros((4,len(ALL)),bool); RES_=np.zeros((4,len(ALL)),bool)
for a,ab in enumerate(ABS):
    s=P[P[f"n_{ab}"].notna()]
    for iso,v in zip(s.isolate.astype(str),s[f"n_{ab}"]>0):
        HAS[a,IX[iso]]=True; RES_[a,IX[iso]]=v
def wlo(k,n,z=1.96):
    if n==0: return 0.0
    p=k/n; dd=1+z*z/n
    return max(0.0,((p+z*z/(2*n))/dd)-z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/dd)
KS=[10,20,30,50,80,120]
PROTECT=None   # optional: Funktion(Merkmalsmenge) -> 0/1 fuer ein vorzeichenbehaftetes Merkmal
def run(ab,use_spec=True,reps=3,dedup=None):
    ai=ABS.index(ab)
    d,f2,bl,gn=build(ab)
    if dedup:
        ks=keep_set(d,META,rule=dedup,ab=ab)
        m=[i for i,iso in enumerate(d.isolate) if iso in ks]
        d=d.iloc[m].reset_index(drop=True); f2=[f2[i] for i in m]
    y=d.R.values; sts=d.st.values; iso=list(d.isolate)
    gidx=np.array([IX[s] for s in iso])
    mcar=collections.defaultdict(list)
    for i,f in enumerate(f2):
        for mm in f: mcar[mm].append(i)
    mcar={mm:np.array(v) for mm,v in mcar.items() if len(v)>=10}
    acc={k:[0,0,0,0] for k in KS}; Pr=np.zeros(len(y))
    for rep in range(reps):
        r=np.random.default_rng(9100+rep); uq=pd.unique(sts).copy(); r.shuffle(uq)
        fo=np.array([{s:i%5 for i,s in enumerate(uq)}[s] for s in sts])
        for f in range(5):
            trm=fo!=f; tr=np.where(trm)[0]; te=np.where(~trm)[0]
            gtr=gidx[tr]; HAStr=HAS[:,gtr]; REStr=RES_[:,gtr]
            pos_in_tr=-np.ones(len(ALL),np.int64); pos_in_tr[gtr]=np.arange(len(gtr))
            tot_n=HAStr.sum(1).astype(float); tot_r=(REStr&HAStr).sum(1).astype(float)
            cand=[]
            for mm,ci in mcar.items():
                sel=ci[trm[ci]]; n=len(sel)
                if n<10: continue
                rr=int(y[sel].sum())
                if use_spec:
                    # Antibiotikaspezifitaet auf der Odds-Ratio-Skala.
                    # Der PPV haengt von der Praevalenz der jeweiligen Substanz ab; ein
                    # Marker mit identischem OR fuer alle vier bekaeme auf der PPV-Skala
                    # 1.6 fuer Ciprofloxacin und 0.6 fuer Tobramycin. Das Odds-Ratio-
                    # Verhaeltnis ist in diesem Fall exakt 1.
                    ci_tr=pos_in_tr[gidx[sel]]; ci_tr=ci_tr[ci_tr>=0]
                    lg=np.empty(4); okall=True
                    for o in range(4):
                        h=HAStr[o,ci_tr]; nn=float(h.sum())
                        if nn<10 or (tot_n[o]-nn)<10: okall=False; break
                        a_=float((REStr[o,ci_tr]&h).sum()); b_=nn-a_
                        c_=tot_r[o]-a_; d_=tot_n[o]-nn-c_
                        lg[o]=np.log(((a_+.5)/(b_+.5))/((c_+.5)/(d_+.5)))
                    if not okall: continue
                    if lg[ai]-(lg.sum()-lg[ai])/3.0 < np.log(1.4): continue
                cand.append((mm,wlo(rr,n),n))
            rk=[mm for mm,_,_ in sorted(cand,key=lambda x:(-x[1],-x[2]))]
            for k in KS:
                Ps=set(rk[:k])
                for i in te:
                    p=bool(f2[i]&Ps); t=y[i]; a=acc[k]
                    if p and t:a[0]+=1
                    elif p:a[1]+=1
                    elif t:a[2]+=1
                    else:a[3]+=1
            T12,T50,T120=set(rk[:12]),set(rk[:50]),set(rk[:120])
            NX=7 if PROTECT is not None else 6
            X=np.zeros((len(f2),NX))
            for i,ff in enumerate(f2):
                X[i,0]=len(ff&T12);X[i,1]=len(ff&T50);X[i,2]=len(ff&T120)
                X[i,3]=any(x.startswith("ACQ:") for x in ff);X[i,4]=any(x.startswith("TRUNC:") for x in ff)
                X[i,5]=any(x.startswith("LOSS:") for x in ff)
                if PROTECT is not None: X[i,6]=PROTECT(ff)
            mdl=LogisticRegression(max_iter=2000).fit(X[tr],y[tr]); Pr[te]+=mdl.predict_proba(X[te])[:,1]/reps
    out={}
    for k in KS:
        TP,FP,FN,TN=acc[k]
        out[k]=dict(sens=TP/max(TP+FN,1),spec=TN/max(TN+FP,1),ppv=TP/max(TP+FP,1),npv=TN/max(TN+FN,1))
    lo=Pr<0.05; hi=Pr>0.80
    br=np.mean((Pr-y)**2); base=np.mean((y.mean()-y)**2)
    return dict(n=len(d),R=int(y.sum()),prev=float(y.mean()),fr=out,bss=float(1-br/base),
        lo=float(lo.mean()),lo_r=float(y[lo].mean()) if lo.sum() else None,
        hi=float(hi.mean()),hi_r=float(y[hi].mean()) if hi.sum() else None,
        Pr=Pr.tolist(), iso=list(d.isolate))
if __name__=="__main__":
    R={}
    for ab in ABS:
        R[ab]={}
        for mode,tag in [(None,"voll"),("first","dedup_first"),("worst","dedup_worst")]:
            x=run(ab,True,3,mode); R[ab][tag]=x
            print(f"{ab:<4}{tag:<14}n={x['n']:>5} R={x['R']:>5} ({x['prev']*100:4.1f}%) "
                  f"k50 {x['fr'][50]['sens']*100:5.1f}/{x['fr'][50]['spec']*100:5.1f} "
                  f"BSS {x['bss']:.3f} lo {x['lo']*100:5.1f}% (R {(x['lo_r'] or 0)*100:.1f}) "
                  f"hi {x['hi']*100:4.1f}% (R {(x['hi_r'] or 0)*100:.1f})",flush=True)
    json.dump(R,open("dedup4.json","w"),indent=1)
