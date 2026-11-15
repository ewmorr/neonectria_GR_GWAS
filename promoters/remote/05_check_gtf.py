#!/usr/bin/env python3
"""Diagnose why StringTie transcripts do not map onto the annotation.
usage: 05_check_gtf.py ref.gff3 sample.gtf [more.gtf ...]
For each transcript (cov >= 1) that overlaps an annotated CDS span, report whether the overlapped gene is on the
SAME or OPPOSITE strand (opposite for most overlaps => the strandedness setting is inverted), and how many
same-strand transcripts contain the gene's whole CDS / its start codon / its stop codon (fragmentation / 3' bias)."""
import sys, re, bisect
from collections import defaultdict
ref, gtfs = sys.argv[1], sys.argv[2:]
by = defaultdict(list)
for l in open(ref):
    if l.startswith("#"): continue
    f = l.rstrip("\n").split("\t")
    if len(f) > 8 and f[2] == "mRNA": by[f[0]].append((int(f[3]), int(f[4]), f[6]))
for c in by: by[c].sort()
maxlen = 20000
for gtf in gtfs:
    n = same = opp = full = five = three = noov = 0
    for l in open(gtf):
        if l.startswith("#"): continue
        f = l.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "transcript": continue
        m = re.search(r'cov "([^"]+)"', f[8])
        if not m or float(m.group(1)) < 1.0: continue
        n += 1
        c, s, e, st = f[0], int(f[3]), int(f[4]), f[6]
        ov = [g for g in by.get(c, []) if g[0] <= e and g[1] >= s] if len(by.get(c, [])) < 5000 else \
             [g for g in by[c][max(0, bisect.bisect_left([x[0] for x in by[c]], s) - 50): bisect.bisect_right([x[0] for x in by[c]], e) + 50] if g[0] <= e and g[1] >= s]
        if not ov: noov += 1; continue
        # strand of the gene with the largest overlap
        g = max(ov, key=lambda g: min(e, g[1]) - max(s, g[0]))
        if g[2] == st:
            same += 1
            if s <= g[0] and e >= g[1]: full += 1
            fiveprime = (s <= g[0]) if st == "+" else (e >= g[1])
            threeprime = (e >= g[1]) if st == "+" else (s <= g[0])
            five += fiveprime; three += threeprime
        else: opp += 1
    ov_t = same + opp
    print(f"{gtf}: transcripts(cov>=1)={n}; no gene overlap={noov}; overlap gene: same strand={same} opposite={opp}"
          + (f" ({100*same/ov_t:.0f}% same)" if ov_t else "")
          + f"; of same-strand: contain full CDS={full}, cover start codon={five}, cover stop codon={three}")
