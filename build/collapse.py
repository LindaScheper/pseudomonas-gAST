# -*- coding: utf-8 -*-
"""Haplotyp-Kollabierung: Marker, die dieselben Isolate treffen, zu einem Block.

Positionen innerhalb eines Locus sind oft vollstaendig gekoppelt — jedes Isolat,
das die eine Variante traegt, traegt auch die andere. Solche Marker sind keine
unabhaengige Information; ein Panel fester Groesse wuerde sich sonst mit
Kopien desselben Haplotyps fuellen. Marker, deren Traegermengen sich mit einem
Jaccard-Koeffizienten von mindestens `thr` ueberlappen, werden deshalb zu einem
Block zusammengefasst.

Regeln
  - nur Marker mit mindestens `minn` Traegern nehmen an der Zusammenfassung teil;
    seltenere bleiben unveraendert stehen
  - zusammengefasst wird ausschliesslich innerhalb eines Locus (Einzelverkettung
    ueber alle Paare des Locus)
  - Namensgeber ist der Marker mit den meisten Traegern, bei Gleichstand der
    alphabetisch erste; der Blockname lautet "<Namensgeber> [+k]" mit k weiteren
    Mitgliedern
  - Traegermenge des Blocks ist die Vereinigung seiner Mitglieder

Rueckgabe: (f2, blocks)
  f2      wie feat, aber mit Blocknamen statt Markernamen
  blocks  Blockname -> Liste der Mitglieder

Geprueft gegen buildcache.pkl des Manuskriptlaufs: fuer alle vier Substanzen
identische Blockmengen und identische Merkmalsmengen je Isolat.
"""
import collections
import itertools

_PREFIX = ("TRUNC:", "LOSS:", "FS:", "STOP:", "ELONG:", "ACQ:")


def _locus(marker):
    for p in _PREFIX:
        if marker.startswith(p):
            return marker[len(p):]
    return marker.rsplit("_", 1)[0]


def collapse_features(feat, thr=0.90, minn=5):
    car = collections.defaultdict(set)
    for i, f in enumerate(feat):
        for m in f:
            car[m].add(i)
    big = {m: s for m, s in car.items() if len(s) >= minn}

    parent = {m: m for m in big}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    by_locus = collections.defaultdict(list)
    for m in big:
        by_locus[_locus(m)].append(m)

    for members in by_locus.values():
        members.sort()
        for a, b in itertools.combinations(members, 2):
            ra, rb = find(a), find(b)
            if ra == rb:
                continue
            A, B = big[a], big[b]
            n = len(A & B)
            if n and n / (len(A) + len(B) - n) >= thr:
                parent[ra] = rb

    comp = collections.defaultdict(list)
    for m in big:
        comp[find(m)].append(m)

    name, blocks = {}, {}
    for members in comp.values():
        rep = min(members, key=lambda m: (-len(big[m]), m))
        nm = rep if len(members) == 1 else "%s [+%d]" % (rep, len(members) - 1)
        blocks[nm] = sorted(members)
        for x in members:
            name[x] = nm

    f2 = [{name.get(m, m) for m in f} for f in feat]
    return f2, blocks


if __name__ == "__main__":
    demo = [{"g_1", "g_2", "g_9"}, {"g_1", "g_2"}, {"g_1", "g_2"},
            {"g_1", "g_2"}, {"g_1", "g_2"}, {"g_9"}]
    f2, bl = collapse_features(demo, 0.90, 5)
    print(f2[0], bl)
