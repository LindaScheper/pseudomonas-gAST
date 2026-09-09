"""Deduplikation: ein Isolat je Patient x Sequenztyp.
Auswahlregel primaer: chronologisch erstes Isolat (CLSI-M39-Logik).
Robustheitsvariante: das resistenteste Isolat der Gruppe (worst case).
Isolate ohne Patientenkennung bleiben einzeln erhalten."""
import pandas as pd, numpy as np

def keep_set(d, meta, rule="first", ab=None):
    pat=[]; st=[]; dt=[]
    for i in d.isolate:
        m=meta.get(i,{})
        p=m.get("patient"); p="__none__%s"%i if p in (None,"","None") or pd.isna(p) else str(p)
        pat.append(p); st.append(str(m.get("mlst","NA")))
        try: dt.append(pd.to_datetime(m.get("date"),errors="coerce"))
        except Exception: dt.append(pd.NaT)
    t=d.copy(); t["_p"]=pat; t["_s"]=st; t["_d"]=dt
    t["_d"]=t["_d"].fillna(pd.Timestamp("2100-01-01"))
    if rule=="first":
        t=t.sort_values(["_p","_s","_d","isolate"])
    else:
        t=t.sort_values(["_p","_s",f"n_{ab}","isolate"],ascending=[True,True,False,True])
    return set(t.drop_duplicates(["_p","_s"]).isolate)
