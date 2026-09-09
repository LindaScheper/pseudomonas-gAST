"""build() fuer alle vier Substanzen einmal rechnen und cachen -> buildcache.pkl"""
import sys, time; sys.path.insert(0, ".")
t0 = time.time()
exec(open("03_core.py").read())
import pickle
print("core geladen", round(time.time() - t0, 1), flush=True)
C = {}
for ab in ["CIP", "CAZ", "MEM", "TOB"]:
    d, f2, bl, gn = build(ab)
    C[ab] = (d.reset_index(drop=True), f2, gn)
    print(ab, len(d), round(time.time() - t0, 1), flush=True)
pickle.dump(C, open("buildcache.pkl", "wb"))
print("buildcache.pkl geschrieben", round(time.time() - t0, 1), flush=True)
pickle.dump(dict(META=META, cov=cov, flip=flip, major=major), open("corebits.pkl", "wb"))
print("fertig", round(time.time() - t0, 1))
