"""Abbildungsdaten fuer die Effluxpumpen-Abbildung."""
import json, pickle, numpy as np, pandas as pd, collections
J={}
J["matrix"]=json.load(open("pumps.json"))
J["rescue"]=json.load(open("porin_rescue.json"))
# Hinweis: J["model"] wird von fig5_pumps.py nicht aus pumpfig.json gelesen,
# sondern direkt aus protect_full.json ueberschrieben (Schluessel voll|ohne,
# voll|mit). Ein eigener Eintrag hier ist daher nicht erforderlich.
# MIC-Modell nachrechnen und ablegen
import sys; sys.path.insert(0,".")
from sklearn.linear_model import Ridge
C=pickle.load(open("buildcache.pkl","rb"))
d,f2,gn=C["MEM"]
A=d.centre.values=="A"; mic=d.v_MEM.values.astype(float); ok=A&np.isfinite(mic)
PORIN=np.array([("TRUNC:oprD" in f) or ("LOSS:oprD" in f) for f in f2])
PROT=np.array([any(m in f for m in ("TRUNC:mexA","LOSS:mexA","TRUNC:mexB","LOSS:mexB")) for f in f2])
P=PROT & ~PORIN
idx=np.where(ok)[0]; sts=d.st.values[idx]; yv=np.log2(mic[idx])
def run(extra):
    reps=3; pred=np.zeros(len(idx))
    for rep in range(reps):
        r=np.random.default_rng(9100+rep); uq=pd.unique(sts).copy(); r.shuffle(uq)
        fo=np.array([{s:i%5 for i,s in enumerate(uq)}[s] for s in sts])
        for f in range(5):
            trm=fo!=f; tr=np.where(trm)[0]; te=np.where(~trm)[0]
            mcar=collections.defaultdict(list)
            for k,i in enumerate(idx):
                for mm in f2[i]: mcar[mm].append(k)
            cand=[]
            for mm,ci in mcar.items():
                ci=np.array(ci); sel=ci[trm[ci]]
                if len(sel)<10: continue
                cand.append((mm, yv[sel].mean(), len(sel)))
            rk=[mm for mm,_,_ in sorted(cand,key=lambda t:(-t[1],-t[2]))]
            T=[set(rk[:12]),set(rk[:50]),set(rk[:120])]
            X=np.zeros((len(idx),6+(1 if extra is not None else 0)))
            for k,i in enumerate(idx):
                ff=f2[i]
                for c in range(3): X[k,c]=len(ff&T[c])
                X[k,3]=any(m.startswith("ACQ:") for m in ff)
                X[k,4]=any(m.startswith("TRUNC:") for m in ff)
                X[k,5]=any(m.startswith("LOSS:") for m in ff)
            if extra is not None: X[:,6]=extra[idx]
            mdl=Ridge(alpha=1.0).fit(X[tr],yv[tr]); pred[te]+=mdl.predict(X[te])/reps
    return dict(rmse=float(np.sqrt(np.mean((pred-yv)**2))), r2=float(1-np.mean((pred-yv)**2)/np.var(yv)))
J["mic"]={"Basismodell":run(None), "mit Marker":run(P)}
J["mic_dist"]={}
lm=mic[ok]; g=P[ok]
for v in sorted(set(lm)):
    m=lm==v
    J["mic_dist"][str(v)]=dict(n=int(m.sum()), n_prot=int((m&g).sum()))
json.dump(J,open("pumpfig.json","w"),indent=1)
print("pumpfig.json geschrieben")
print(" MIC-Modell:", {k:(round(v['r2'],3),round(v['rmse'],2)) for k,v in J["mic"].items()})
