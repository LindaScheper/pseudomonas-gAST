"""Pfade an einer Stelle.

DATA  Verzeichnis mit den vorberechneten Eingaben (figdata.pkl, specratio.pkl,
      specratio.json, fs_fig.json, final_all.json, elong.json)
OUT   Verzeichnis, in das die Abbildungen geschrieben werden
"""
import os
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.environ.get("GAST_DATA", HERE) + os.sep
OUT  = os.environ.get("GAST_OUT",  "/home/sch22/analysis_scripts/analysis") + os.sep
os.makedirs(OUT, exist_ok=True)