"""Kennzahlen fuer die sieben Publikationsabbildungen -> figdata.pkl (nutzt buildcache.pkl)"""
import numpy as np, pandas as pd, collections, pickle, json, re, time
t0 = time.time()
C = pickle.load(open("buildcache.pkl", "rb"))
CB = pickle.load(open("corebits.pkl", "rb"))
META, cov, flip, major = CB["META"], CB["cov"], CB["flip"], CB["major"]
ABS = ["CIP", "CAZ", "MEM", "TOB"]
D = {}

def det_set(f2, y, min_n=10, min_ppv=0.5):
    car = collections.defaultdict(list)
    for i, f in enumerate(f2):
        for m in f: car[m].append(i)
    DET = {m for m, c in car.items() if len(c) >= min_n and y[c].mean() >= min_ppv}
    DET |= {m for m in car if m.startswith("ACQ:")}
    return DET, car

# ---------- 1 Kohorte
coh = {}
for ab in ABS:
    d, f2, gn = C[ab]
    coh[ab] = dict(n=len(d), R=int(d.R.sum()),
                   norm_A=d.loc[d.centre == "A", "norm"].values,
                   norm_B=d.loc[d.centre == "B", "norm"].values,
                   nA=int((d.centre == "A").sum()), nB=int((d.centre == "B").sum()),
                   rA=float(d.loc[d.centre == "A", "R"].mean()),
                   rB=float(d.loc[d.centre == "B", "R"].mean()))
D["cohort"] = coh
import pandas as pd, numpy as np, collections

# d0 bleibt fuer die Jahres-Uebersicht (D["years"]) bei CIP, da dort das
# vollstaendige Datumsfeld am zuverlaessigsten vorliegt; falls gewuenscht,
# kann das ebenfalls auf die Union umgestellt werden.
d0 = C["CIP"][0]
yrs = []
for i in d0.isolate:
    try: yrs.append(pd.to_datetime(META.get(i, {}).get("date")).year)
    except Exception: yrs.append(np.nan)
tmp = d0.assign(year=yrs)
D["years"] = tmp.dropna(subset=["year"]).groupby(["year", "centre"]).size().unstack(fill_value=0)

# --- Sequenztyp-Verteilung: Union ueber ALLE Antibiotika-Kohorten (5749 Isolate),
# nicht nur CIP. Jedes Isolat hat unabhaengig vom Antibiotikum denselben ST,
# ein einfaches Update-Mapping reicht daher aus.
_st_map = {}
for _ab in ABS:
    _d, _f2, _gn = C[_ab]
    for _iso, _st in zip(_d.isolate, _d.st):
        _st_map[_iso] = _st

_st_counts = collections.Counter(_st_map.values())
D["st_top"] = _st_counts.most_common(16)
D["st_all"] = dict(_st_counts)
D["n_cohort_union"] = len(_st_map)   # 5749, zur Kontrolle/Beschriftung nutzbar
D["highrisk"] = {"235","111","175","244","277","357","308","381","233","654"}
D["npat"] = {ab: int(pd.Series([str(META.get(i, {}).get("patient")) for i in C[ab][0].isolate]).nunique())
             for ab in ABS}
print("1 ok", round(time.time()-t0,1), flush=True)

# ---------- 2 Referenzpolarisierung
M = pd.read_pickle("v3_mut.pkl")
M = M[~M.gene.isin({"mexA","mexB","mexY","mexT","colR"})]
POS = re.compile(r"(\d+)")
alt = collections.defaultdict(collections.Counter)
for g, s in zip(M.gene.values, M.mut_str.values):
    if not s: continue
    for m in str(s).split(","):
        m = m.strip()
        if not m: continue
        dd = POS.search(m)
        if dd: alt[(g, int(dd.group(1)))][m] += 1
rows = []
for (g, p), c in alt.items():
    m, k = c.most_common(1)[0]
    rows.append((g, p, m, k / cov[g]))   # Frequenz der haeufigsten Variante an dieser Position
D["varfreq"] = np.array([r[3] for r in rows])
D["nearfixed"] = sorted([r for r in rows if r[3] >= 0.90], key=lambda x: -x[3])[:12]
D["nflip"] = len(flip); D["nmaj"] = len(major)
D["spain_filter"] = [("pmrA",208,6),("ampD",196,6),("parR",122,5),("pmrB",275,17),
                     ("mexX",298,22),("parS",295,27),("mexR",119,16),("armZ",295,44),
                     ("mexZ",100,51),("gyrA",113,67),("oprD",38,35)]
print("2 ok", round(time.time()-t0,1), flush=True)

# ---------- 3 Trunkierung / Genverlust
tr = []
for ab in ABS:
    d, f2, gn = C[ab]; y = d.R.values
    car = collections.defaultdict(list)
    for i, f in enumerate(f2):
        for m in f:
            if m.startswith(("TRUNC:", "LOSS:")): car[m].append(i)
    tot = int(y.sum()); N = len(y)
    for m, c in car.items():
        if len(c) < 15: continue
        n = len(c); r = int(y[c].sum())
        a, b, cc, dd2 = r, n - r, tot - r, (N - n) - (tot - r)
        if min(a, b, cc, dd2) == 0: continue
        orr = (a / b) / (cc / dd2); se = np.sqrt(1/a + 1/b + 1/cc + 1/dd2)
        tr.append(dict(m=m, ab=ab, n=n, ppv=r/n, rest=cc/(cc+dd2), orr=orr,
                       lo=float(np.exp(np.log(orr)-1.96*se)), hi=float(np.exp(np.log(orr)+1.96*se))))
D["trunc"] = tr
print("3 ok", round(time.time()-t0,1), flush=True)

# ---------- 4 / 5b MIC-stratifizierte Detektion
det = {}
for ab in ABS:
    d, f2, gn = C[ab]; y = d.R.values
    DET, car = det_set(f2, y)
    has = np.array([bool(f & DET) for f in f2])
    A = d[d.type == "mic_V"]
    det[ab] = dict(mic=[(float(v), len(g), float(has[g.index].mean()),
                         float((g[f"n_{ab}"] > 0).mean()))
                        for v, g in A.groupby(f"v_{ab}") if len(g) >= 20],
                   bp=float(A[f"bp_{ab}"].iloc[0]),
                   sens_overall=float(has[y].mean()),
                   floor=float(has[(d.type == "mic_V") & (~d.R)].mean()))
    if ab == "CIP":
        cn = {m for m in car if m.rsplit("_",1)[0] in ("gyrA","parC")
              and m.rsplit("_",1)[1].split(" ")[0] in ("83","87","80")}
        hc = np.array([bool(f & cn) for f in f2])
        det[ab]["canon"] = [(float(v), len(g), float(hc[g.index].mean()))
                            for v, g in A.groupby(f"v_{ab}") if len(g) >= 20]
        det[ab]["canon_overall"] = float(hc[y].mean())
D["detect"] = det
print("4 ok", round(time.time()-t0,1), flush=True)

# ---------- 5 Frontier
_fa = json.load(open("final_all.json"))
D["frontier"] = {ab: {"mit": _fa[ab]["voll"]} for ab in ABS}
D["specratio"] = [("PBP4 trunc. → CAZ",5.46),("ampD trunc. → CAZ",3.32),("fusA1 671 → TOB",3.00),
                  ("gyrA S83 → CIP",2.23),("parC S87 → CIP",1.50),("oprD trunc. → MEM",1.39),
                  ("sltB2 293",1.22),("galU 37",1.11),("oprJ 25",1.09),
                  ("ampG 439",0.95),("mexZ 126",0.93),("amgK 154",0.89)]

# ---------- 5b CAZ-Mechanismusklassen
d, f2, gn = C["CAZ"]; y = d.R.values; R = int(y.sum())
DET, car = det_set(f2, y)
AMPC = {"ampC","ampD","ampDh2","ampDh3","ampR","ampG","ampO","ampP","mpl","nagZ","nagA","anmK","dacC","creB","creC","creD"}
PBPS = {g for g in gn if g.startswith("PBP")}
EFF = {"mexR","nalC","nalD","oprM"}
def gene_of(m):
    if m.startswith("ACQ:"): return "ACQ"
    if m.startswith(("TRUNC:","LOSS:")): return m.split(":")[1]
    return m.rsplit("_",1)[0]
cls = []
for lab, S in [("Acquired enzymes",{"ACQ"}),("AmpC pathway",AMPC),("PBPs",PBPS),
               ("Efflux regulation",EFF),("Porin / LPS",{"oprD","galU"}),
               ("Cell wall, other",set(gn)-AMPC-PBPS-EFF-{"oprD","galU"})]:
    ms = [m for m in DET if gene_of(m) in S]
    cc = set()
    for m in ms: cc |= set(car[m])
    cc = sorted(cc)
    if cc: cls.append(dict(lab=lab, n=len(ms), carr=len(cc),
                           ppv=float(y[cc].mean()), cover=float(y[cc].sum()/R)))
D["caz_classes"] = cls
print("5 ok", round(time.time()-t0,1), flush=True)

# ---------- 6 Deduplikation und 7 Befundzustaende: direkt aus dem Endlauf
FA = json.load(open("final_all.json"))
D["final"] = FA
D["dedup"] = {}
for ab in ABS:
    v = FA[ab]
    D["dedup"][ab] = dict(
        prev=[v[t]["prev"]*100 for t in ("voll","dedup_first","dedup_worst")],
        bss=[v[t]["bss"] for t in ("voll","dedup_first","dedup_worst")],
        sens=[v["voll"]["fr"]["50"]["sens"]*100, v["dedup_first"]["fr"]["50"]["sens"]*100],
        spec=[v["voll"]["fr"]["50"]["spec"]*100, v["dedup_first"]["fr"]["50"]["spec"]*100])
D["report"] = {ab: (FA[ab]["voll"]["lo"]*100,
                    (FA[ab]["voll"]["lo_r"] or 0)*100 if FA[ab]["voll"]["lo"] > 0 else None,
                    FA[ab]["voll"]["hi"]*100,
                    (FA[ab]["voll"]["hi_r"] or 0)*100) for ab in ABS}
D["bss"] = {ab: FA[ab]["voll"]["bss"] for ab in ABS}
D["ngenes"] = {ab: len(C[ab][2]) for ab in ABS}
D["dedup_markers"] = [("parC S87",99.6,98.7,76),("gyrA S83",94.4,88.2,211),
                      ("rmtB4",100,100,6),("aac(6')-Ib",93.8,100,5),
                      ("ampD trunc.",67.7,65.2,23),("PBP4 trunc.",85.7,71.4,14),
                      ("oprD trunc.",48.0,44.0,166),("nalD trunc.",67.4,45.2,31),
                      ("mexZ pos.126",84.8,63.6,11),("fusA1 pos.154",78.8,33.3,9)]

pickle.dump(D, open("figdata.pkl", "wb"))
print("figdata.pkl geschrieben", round(time.time()-t0,1))
for k, v in D.items():
    print(f"  {k:<16}{type(v).__name__} {len(v) if hasattr(v,'__len__') else ''}")