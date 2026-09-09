"""Konfidenzintervalle auf den Kernzahlen des Manuskripts.

Der Brier Skill Score wird mit einem Cluster-Bootstrap ueber Sequenztypen
intervalliert: gezogen werden ganze Sequenztypen, nicht einzelne Isolate, weil
Isolate einer Linie nicht unabhaengig sind. Die gepaarte Differenz zwischen zwei
Modellen wird auf denselben Ziehungen gebildet und ist dadurch deutlich enger
als die Differenz der beiden Randintervalle.

Anteile (Restresistenz, Rule-in-Praezision) bekommen Wilson-Intervalle.

Erzeugt daten/ci_kern.json und schreibt die Tabelle nach stdout.
"""
import sys, os, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "vorstufen"))
sys.path.insert(0, os.path.join(HERE, ".."))
from paths import DATA

ABS = ["CIP", "CAZ", "MEM", "TOB"]
ON  = json.load(open(DATA + "final_all.json"))
OFF = json.load(open(DATA + "spec_off.json"))
rng = np.random.default_rng(7)


def wilson(k, n, z=1.96):
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (100 * (c - h), 100 * (c + h))


def bss(y, p):
    return 1 - np.mean((p - y) ** 2) / np.mean((y.mean() - y) ** 2)


def clusters(stv):
    g = {}
    for i, s in enumerate(stv):
        g.setdefault(s, []).append(i)
    return g


def boot(y, stv, fn, reps=2000):
    g = clusters(stv); keys = list(g); out = []
    for _ in range(reps):
        sel = rng.choice(len(keys), len(keys), replace=True)
        idx = np.concatenate([g[keys[k]] for k in sel])
        if y[idx].min() == y[idx].max():
            continue
        out.append(fn(idx))
    return np.percentile(out, 2.5), np.percentile(out, 97.5), np.array(out)


import pickle
C = pickle.load(open(DATA + "buildcache.pkl", "rb"))
OUT = {}
for ab in ABS:
    d = C[ab][0]
    ymap = dict(zip(d.isolate, d.R.values.astype(float)))
    stm  = dict(zip(d.isolate, d.st.values))
    A, B = ON[ab]["voll"], OFF[ab]
    iso = A["iso"]
    y   = np.array([ymap[i] for i in iso], float)
    stv = np.array([stm[i] for i in iso])
    pa, pb = np.array(A["Pr"]), np.array(B["Pr"])

    lo_a, hi_a, _ = boot(y, stv, lambda ix: bss(y[ix], pa[ix]))
    lo_d, hi_d, ds = boot(y, stv, lambda ix: bss(y[ix], pa[ix]) - bss(y[ix], pb[ix]))
    lo, hi = pa < 0.05, pa > 0.80
    OUT[ab] = dict(
        bss=bss(y, pa), bss_lo=lo_a, bss_hi=hi_a,
        dbss=bss(y, pa) - bss(y, pb), dbss_lo=lo_d, dbss_hi=hi_d,
        p_gt0=float((ds > 0).mean()),
        S_share=100 * lo.mean(), S_n=int(lo.sum()),
        S_r=100 * y[lo].mean() if lo.sum() else None, S_ci=wilson(y[lo].sum(), lo.sum()),
        R_share=100 * hi.mean(), R_n=int(hi.sum()),
        R_r=100 * y[hi].mean() if hi.sum() else None, R_ci=wilson(y[hi].sum(), hi.sum()))
    o = OUT[ab]
    print("%-4s BSS %.3f [%.3f-%.3f]   Filtergewinn %+.3f [%+.3f bis %+.3f], P(>0) %.2f" % (
        ab, o["bss"], o["bss_lo"], o["bss_hi"], o["dbss"], o["dbss_lo"], o["dbss_hi"], o["p_gt0"]))
    print("     nicht-resistent %.1f %% (n=%d), darin R %.1f %% [%.1f-%.1f]" % (
        o["S_share"], o["S_n"], o["S_r"] or 0, o["S_ci"][0], o["S_ci"][1]))
    print("     resistent       %.1f %% (n=%d), darin R %.1f %% [%.1f-%.1f]" % (
        o["R_share"], o["R_n"], o["R_r"] or 0, o["R_ci"][0], o["R_ci"][1]))

# --- Das vorzeichenbehaftete Pumpenmerkmal bei Meropenem -------------------
if os.path.exists(DATA + "noprot.json"):
    A = ON["MEM"]["voll"]
    Bn = json.load(open(DATA + "noprot.json"))["MEM"]
    d = C["MEM"][0]
    ymap = dict(zip(d.isolate, d.R.values.astype(float)))
    stm  = dict(zip(d.isolate, d.st.values))
    iso = A["iso"]
    assert iso == Bn["iso"], "unterschiedliche Isolatreihenfolge"
    y   = np.array([ymap[i] for i in iso], float)
    stv = np.array([stm[i] for i in iso])
    pa, pb = np.array(A["Pr"]), np.array(Bn["Pr"])
    lo_p, hi_p, dp = boot(y, stv, lambda ix: bss(y[ix], pa[ix]) - bss(y[ix], pb[ix]))
    OUT["MEM"]["pump_dbss"] = bss(y, pa) - bss(y, pb)
    OUT["MEM"]["pump_lo"], OUT["MEM"]["pump_hi"] = lo_p, hi_p
    OUT["MEM"]["pump_p_gt0"] = float((dp > 0).mean())
    print("\nPumpenmerkmal, Meropenem: dBSS %+.3f [%+.3f bis %+.3f], P(>0) %.2f" % (
        OUT["MEM"]["pump_dbss"], lo_p, hi_p, OUT["MEM"]["pump_p_gt0"]))
    for tag, p in [("mit ", pa), ("ohne", pb)]:
        s = p < 0.05
        print("   %s nicht-resistent %.1f %% (n=%d), darin R %.1f %%" % (
            tag, 100 * s.mean(), s.sum(), 100 * y[s].mean()))

# --- Das MexXY-Merkmal bei Tobramycin ---------------------------------------
if os.path.exists(DATA + "tob_xyprot.json"):
    A = json.load(open(DATA + "tob_xyprot.json"))["TOB"]
    Bn = ON["TOB"]["voll"]
    d = C["TOB"][0]
    ymap = dict(zip(d.isolate, d.R.values.astype(float)))
    stm  = dict(zip(d.isolate, d.st.values))
    iso = A["iso"]
    assert iso == Bn["iso"], "unterschiedliche Isolatreihenfolge"
    y   = np.array([ymap[i] for i in iso], float)
    stv = np.array([stm[i] for i in iso])
    pa, pb = np.array(A["Pr"]), np.array(Bn["Pr"])
    lo_x, hi_x, dx = boot(y, stv, lambda ix: bss(y[ix], pa[ix]) - bss(y[ix], pb[ix]))
    OUT["TOB"]["mexxy_dbss"] = bss(y, pa) - bss(y, pb)
    OUT["TOB"]["mexxy_lo"], OUT["TOB"]["mexxy_hi"] = lo_x, hi_x
    OUT["TOB"]["mexxy_p_gt0"] = float((dx > 0).mean())
    print("\nMexXY-Merkmal, Tobramycin: BSS %.3f -> %.3f, dBSS %+.3f [%+.3f bis %+.3f], P(>0) %.2f" % (
        bss(y, pb), bss(y, pa), OUT["TOB"]["mexxy_dbss"], lo_x, hi_x, OUT["TOB"]["mexxy_p_gt0"]))

json.dump(OUT, open(DATA + "ci_kern.json", "w"), indent=1, default=float)
print("\ngeschrieben:", DATA + "ci_kern.json")
