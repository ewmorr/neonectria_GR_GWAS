#!/usr/bin/env python3
"""Infer library strandedness from SAM on stdin (proper pairs, primary), using non-overlapping genes from a gff3.
usage: samtools view -f 2 -F 256 x.bam | check_strand.py genes.gff3
"""
import sys, bisect, re
from collections import defaultdict

genes = defaultdict(list)
for line in open(sys.argv[1]):
    if line.startswith("#"): continue
    f = line.rstrip("\n").split("\t")
    if len(f) > 8 and f[2] == "mRNA":
        genes[f[0]].append((int(f[3]), int(f[4]), f[6]))
idx = {}
for c, g in genes.items():
    g.sort()
    keep = []
    run_max = -1                              # running max end of all previous genes
    for i, (s, e, st) in enumerate(g):       # drop genes overlapping a neighbour (strand ambiguity)
        nxt_start = g[i+1][0] if i + 1 < len(g) else 10**12
        if run_max < s and nxt_start > e: keep.append((s, e, st))
        run_max = max(run_max, e)
    idx[c] = ([k[0] for k in keep], keep)

same = opp = n_total = 0
for line in sys.stdin:
    f = line.split("\t")
    flag = int(f[1])
    if not flag & 0x40: continue                       # read1 only
    if "NH:i:1" not in line: continue                  # unique
    n_total += 1
    c, pos = f[2], int(f[3])
    if c not in idx: continue
    starts, keep = idx[c]
    i = bisect.bisect_right(starts, pos) - 1
    if i < 0 or not (keep[i][0] <= pos <= keep[i][1]): continue
    rstrand = "-" if flag & 0x10 else "+"
    if rstrand == keep[i][2]: same += 1
    else: opp += 1

tot = same + opp
if tot < 1000:
    sys.exit(f"too few informative reads ({tot}); check mapping")
fs = same / tot
print(f"read1 informative: {tot} of {n_total} unique read1; same strand as gene: {fs:.3f}; opposite: {1-fs:.3f}")
if fs > 0.85:   print("=> STRANDED=FR   (read1 sense; fr-secondstrand)")
elif fs < 0.15: print("=> STRANDED=RF   (read1 antisense; dUTP / fr-firststrand)")
else:           print("=> STRANDED=none (unstranded)")
