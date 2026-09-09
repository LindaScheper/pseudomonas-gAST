"""Befundschema mit einem markerbasierten Rule-in-Pfad.

Die Rule-in-Liste ist vorab aus der Substratspezifitaet der Enzymklassen
festgelegt, nicht aus den eigenen praediktiven Werten abgeleitet:
  Meropenem   Metallo-Beta-Laktamasen (IMP, VIM, NDM, FIM, SPM, GIM) und
              Klasse-A-Carbapenemasen (GES-5, GES-13, KPC)
  Ceftazidim  dieselben, dazu die Extended-Spectrum-Enzyme (PER, VEB, PME,
              BEL, GES-1, OXA-10-Gruppe mit erweitertem Spektrum: OXA-10,
              OXA-14/-15/-17/-35/-141)
  Tobramycin  16S-rRNA-Methyltransferasen (ArmA, RmtB) und die
              Aminoglykosid-modifizierenden Enzyme mit Tobramycin im
              Substratspektrum: AAC(6')-I, ANT(4')-II, ANT(2'')-I
"""
import pickle, numpy as np, collections, json, sys; sys.path.insert(0,".")
import dedup_run as D
C=pickle.load(open("buildcache.pkl","rb"))
D.build=lambda ab: (C[ab][0], C[ab][1], {}, C[ab][2])

MBL   = ("blaIMP","blaVIM","blaNDM","blaFIM","blaSPM","blaGIM")
CARBA = MBL + ("blaGES-5","blaGES-13","blaKPC")
ESBL  = ("blaPER","blaVEB","blaPME","blaBEL","blaGES-1","blaOXA-10","blaOXA-14",
         "blaOXA-15","blaOXA-17","blaOXA-35","blaOXA-141","blaSHV-2")
RULEIN={"MEM": CARBA, "CAZ": CARBA+ESBL,
        "TOB": ("armA","rmtB","aac(6')-I","ant(4')-II","ant(2'')-I"), "CIP": ()}
def is_rulein(m, ab):
    if not m.startswith("ACQ:"): return False
    g=m.replace("ACQ:","")
    return any(g.startswith(p) for p in RULEIN[ab])

OUT={}
print(f"{'':<5}{'Rule-in-Traeger':>16}{'% Kohorte':>11}{'davon R':>9}{'Anteil aller R':>16}"
      f"{'davon vom Modell blank':>24}")
print("-"*82)
for ab in ["CIP","CAZ","MEM","TOB"]:
    x=D.run(ab,True,3,None)
    d,f2,gn=C[ab]; y=d.R.values; Pr=np.array(x["Pr"])
    ri=np.array([any(is_rulein(m,ab) for m in f) for f in f2])
    lo=Pr<0.05; hi=Pr>0.80
    blank_before = (~lo)&(~hi)
    if ri.sum():
        print(f"{ab:<5}{int(ri.sum()):>16}{ri.mean()*100:>10.1f}%{y[ri].mean()*100:>8.1f}%"
              f"{y[ri].sum()/y.sum()*100:>15.1f}%{int((ri&blank_before).sum()):>24}")
    else:
        print(f"{ab:<5}{'—':>16}   (keine erworbenen Gene mit Rule-in-Status)")
    # neues Schema
    R2 = hi | ri
    S2 = lo & ~ri
    B2 = ~R2 & ~S2
    OUT[ab]=dict(
        ri_n=int(ri.sum()), ri_cov=float(ri.mean()), ri_ppv=float(y[ri].mean()) if ri.sum() else None,
        ri_share=float(y[ri].sum()/y.sum()) if ri.sum() else 0.0,
        ri_rescued=int((ri&blank_before).sum()),
        old=dict(S=float(lo.mean()), S_r=float(y[lo].mean()) if lo.sum() else None,
                 R=float(hi.mean()), R_r=float(y[hi].mean()) if hi.sum() else None,
                 B=float(blank_before.mean())),
        new=dict(S=float(S2.mean()), S_r=float(y[S2].mean()) if S2.sum() else None,
                 R=float(R2.mean()), R_r=float(y[R2].mean()) if R2.sum() else None,
                 B=float(B2.mean())))
json.dump(OUT,open("rulein_final.json","w"),indent=1)
print("\n=== Befundzustaende: bisher (nur Modell) gegen neu (Modell + Enzym-Rule-in) ===")
print(f"{'':<5}{'sensibel':>20}{'resistent':>22}{'blank':>16}")
for ab in ["CIP","CAZ","MEM","TOB"]:
    o=OUT[ab]["old"]; n=OUT[ab]["new"]
    f=lambda v: "—" if v is None else f"{v*100:.1f}"
    print(f"{ab:<5}{f(o['S']):>8}->{f(n['S']):<8}   {f(o['R'])}(R {f(o['R_r'])}) -> {f(n['R'])}(R {f(n['R_r'])})"
          f"   {f(o['B'])} -> {f(n['B'])}")
