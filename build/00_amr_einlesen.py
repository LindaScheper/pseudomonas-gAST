# -*- coding: utf-8 -*-
"""AMRFinderPlus-Treffer einlesen  ->  v3_amr.pkl

Die Spalten werden auf die kurzen Namen abgebildet, die im uebrigen Code
verwendet werden:

    Isolat                    -> isolate
    Element symbol            -> gene
    Scope / Type / Subtype    -> scope / typ / subtyp
    Class / Subclass          -> cls / subcls
    % Coverage of reference   -> cov
    % Identity to reference   -> ident

Alle uebrigen Spalten bleiben unveraendert erhalten, damit spaetere Filter
(Method, Alignment length, HMM accession) weiter moeglich sind.

Eingabe : daten/AMRFP_all_isolates.xlsx
Ausgabe : daten/v3_amr.pkl
"""
import os, sys
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
from paths import DATA

XLS = DATA + "AMRFP_all_isolates.xlsx"
REN = {"Isolat": "isolate", "Element symbol": "gene", "Scope": "scope",
       "Type": "typ", "Subtype": "subtyp", "Class": "cls", "Subclass": "subcls",
       "% Coverage of reference": "cov", "% Identity to reference": "ident"}

A = pd.read_excel(XLS, sheet_name=0)
A = A.rename(columns=REN)
A["isolate"] = A.isolate.astype(str).str.strip()
A["gene"] = A.gene.astype(str).str.strip()
A.to_pickle(DATA + "v3_amr.pkl")

print("v3_amr.pkl: %d Treffer, %d Isolate, %d verschiedene Genbezeichner"
      % (len(A), A.isolate.nunique(), A.gene.nunique()))
print("Typen:", dict(A.typ.value_counts()))
print("nur erworbene Gene (scope == 'plus' oder typ == 'AMR' ohne EFFLUX): %d"
      % int(((A.typ == "AMR") & (A.cls != "EFFLUX")).sum()))
print("geschrieben:", DATA + "v3_amr.pkl")