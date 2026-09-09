"""Proteinverlaengerung als bisher nicht gewertete Markerklasse.

Trunkierung wird gewertet (rel < 0.95), Verlaengerung (rel > 1.05) bisher nicht.
Eine Insertion oder ein Frameshift, der ueber das Stopcodon hinausliest, kann
ein Protein genauso inaktivieren wie ein vorzeitiges Stopcodon.
"""
import pandas as pd, numpy as np, pickle, collections, json
from scipy.stats import fisher_exact
import sys; sys.path.insert(0, ".")
from dedup import keep_set

M = pd.read_pickle("v3_mut.pkl"); M = M[~M.gene.isin({"mexT", "colR"})]
MODE = M.dropna(subset=["plen"]).groupby("gene").plen.agg(lambda s: s.mode().iloc[0])
M = M.dropna(subset=["plen"]).copy(); M["rel"] = M.plen / M.gene.map(MODE)
LONG = M[M.rel > 1.05]
carr = {g: set(s.isolate) for g, s in LONG.groupby("gene")}

C = pickle.load(open("buildcache.pkl", "rb"))
CB = pickle.load(open("corebits.pkl", "rb")); META = CB["META"]
OUT = {}

print(f"Verlaengerungsaufrufe: {len(LONG):,} in {LONG.isolate.nunique():,} Isolaten "
      f"(Trunkierung: {(M.rel<0.95).sum():,} in {M[M.rel<0.95].isolate.nunique():,})")
print(f"\n{'Gen':<10}{'AB':<5}{'n':>6}{'%R':>8}{'Rest':>8}{'OR':>8}{'95%-KI':>16}{'p':>10}{'MH(ST)':>9}{'dedup n':>9}{'dedup %R':>10}")
rows = []
for ab in ["CIP", "CAZ", "MEM", "TOB"]:
    d, f2, gn = C[ab]; y = d.R.values
    iso2i = {iso: i for i, iso in enumerate(d.isolate)}
    ks = keep_set(d, META, rule="first")
    for g, S in carr.items():
        ii = [iso2i[x] for x in S if x in iso2i]
        if len(ii) < 20: continue
        a = int(y[ii].sum()); b = len(ii) - a
        rest = [i for i in range(len(d)) if i not in set(ii)]
        c = int(y[rest].sum()); dq = len(rest) - c
        if min(a, b, c, dq) == 0: continue
        orr = (a/b)/(c/dq); se = np.sqrt(1/a + 1/b + 1/c + 1/dq)
        lo, hi = np.exp(np.log(orr)-1.96*se), np.exp(np.log(orr)+1.96*se)
        _, p = fisher_exact([[a, b], [c, dq]])
        if orr < 2 or p > 1e-4: continue
        # ST-geschichtet
        st = d.st.values; num = den = 0.0; iis = set(ii)
        for s in pd.unique(st):
            jj = [i for i in range(len(d)) if st[i] == s]
            A = sum(1 for i in jj if i in iis and y[i]); B = sum(1 for i in jj if i in iis and not y[i])
            Cc = sum(1 for i in jj if i not in iis and y[i]); Dd = sum(1 for i in jj if i not in iis and not y[i])
            n = A+B+Cc+Dd
            if n and (A+B) and (Cc+Dd): num += A*Dd/n; den += B*Cc/n
        mh = num/den if den > 0 else np.nan
        dd_i = [i for i in ii if d.isolate.iloc[i] in ks]
        dd_r = [i for i in range(len(d)) if d.isolate.iloc[i] in ks and i not in iis]
        print(f"{g:<10}{ab:<5}{len(ii):>6}{a/len(ii)*100:>7.1f}%{c/(c+dq)*100:>7.1f}%{orr:>8.2f}"
              f"{f'{lo:.2f}-{hi:.2f}':>16}{p:>10.1e}{mh:>9.2f}{len(dd_i):>9}"
              f"{(y[dd_i].mean()*100 if dd_i else 0):>9.1f}%")
        rows.append(dict(gene=g, ab=ab, n=len(ii), ppv=a/len(ii), rest=c/(c+dq), orr=orr,
                         lo=lo, hi=hi, p=p, mh=float(mh), n_dd=len(dd_i),
                         ppv_dd=float(y[dd_i].mean()) if dd_i else None,
                         rest_dd=float(y[dd_r].mean())))
OUT["rows"] = rows

print("\n--- Ueberschneidung mit der Trunkierungsklasse (Meropenem, oprD)")
d, f2, gn = C["MEM"]; y = d.R.values
iso2i = {iso: i for i, iso in enumerate(d.isolate)}
lo_ = [iso2i[x] for x in carr.get("oprD", set()) if x in iso2i]
tr = [i for i, f in enumerate(f2) if "TRUNC:oprD" in f]
both = set(lo_) & set(tr)
print(f"  verlaengert {len(lo_)}, trunkiert {len(tr)}, beides {len(both)}")
only_long = [i for i in lo_ if i not in set(tr)]
print(f"  nur verlaengert: n={len(only_long)}, R {y[only_long].mean()*100:.1f}% (Kohorte {y.mean()*100:.1f}%)")

print("\n--- Wie viele zusaetzliche Resistente wuerden erkannt?")
for ab in ["CIP", "CAZ", "MEM", "TOB"]:
    d, f2, gn = C[ab]; y = d.R.values
    iso2i = {iso: i for i, iso in enumerate(d.isolate)}
    car = collections.defaultdict(list)
    for i, f in enumerate(f2):
        for m in f: car[m].append(i)
    DET = {m for m, cc in car.items() if len(cc) >= 10 and y[cc].mean() >= 0.5} | {m for m in car if m.startswith("ACQ:")}
    cov = set()
    for m in DET: cov |= set(car[m])
    add = set()
    for g in ("oprD", "nfxB"):
        ii = [iso2i[x] for x in carr.get(g, set()) if x in iso2i]
        if len(ii) >= 20 and y[ii].mean() >= 0.5: add |= set(ii)
    newR = {i for i in add if y[i] and i not in cov}
    print(f"  {ab}: bisher {int(y[sorted(cov)].sum())}/{int(y.sum())} erkannt "
          f"(+{len(newR)} durch Verlaengerung = +{len(newR)/int(y.sum())*100:.1f} Punkte)")
    OUT.setdefault("gain", {})[ab] = dict(before=int(y[sorted(cov)].sum()), add=len(newR), tot=int(y.sum()))

json.dump(OUT, open("elong.json", "w"), indent=1, default=float)
print("\nelong.json geschrieben")
