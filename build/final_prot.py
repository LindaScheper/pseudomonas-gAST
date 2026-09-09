"""Endlauf mit dem vorzeichenbehafteten Merkmal fuer alle vier Substanzen und drei Modi."""
import sys, json, pickle, os; sys.path.insert(0,".")
import dedup_run as D
C=pickle.load(open("buildcache.pkl","rb"))
D.build=lambda ab: (C[ab][0], C[ab][1], {}, C[ab][2])
PORIN=("TRUNC:oprD","LOSS:oprD"); PUMP=("TRUNC:mexA","LOSS:mexA","TRUNC:mexB","LOSS:mexB")
D.PROTECT=lambda ff: 1.0 if (any(m in ff for m in PUMP) and not any(m in ff for m in PORIN)) else 0.0
F="final_all.json"; R=json.load(open(F))
todo=[(ab,mode,tag) for ab in ["CIP","CAZ","MEM","TOB"]
      for mode,tag in [(None,"voll"),("first","dedup_first"),("worst","dedup_worst")]
      if not R.get(ab,{}).get(tag,{}).get("_prot")]
if not todo: print("nichts offen"); raise SystemExit
ab,mode,tag=todo[0]
x=D.run(ab,True,3,mode); x["_prot"]=True
R.setdefault(ab,{})[tag]=x
json.dump(R,open(F,"w"),indent=1)
print(f"{ab} {tag:<13}BSS {x['bss']:6.3f}  k50 {x['fr'][50]['sens']*100:5.1f}/{x['fr'][50]['spec']*100:5.1f}  "
      f"sensibel {x['lo']*100:5.1f}% (R {(x['lo_r'] or 0)*100:4.1f})  resistent {x['hi']*100:4.1f}% "
      f"(R {(x['hi_r'] or 0)*100:5.1f})   [{len(todo)-1} offen]")
