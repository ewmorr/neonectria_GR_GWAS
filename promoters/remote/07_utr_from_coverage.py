#!/usr/bin/env python3
"""Call 5'/3' UTR extensions from strand-specific, pooled RNA-seq coverage.
usage: 07_utr_from_coverage.py ref.gff3 bamdir out.tsv sample1 sample2 ...   (needs bamdir/<sample>.{plus,minus}.bam, samtools on PATH)

For each gene the window runs outward from the start (5') or stop (3') codon, clipped at the nearest neighbouring CDS
(any strand) and at MAX_WIN. Reference expression = median depth of the 150 bp of CDS next to that edge. The UTR end
is the furthest point reachable from the CDS edge with 21-bp-smoothed depth >= max(MIN_ABS, FRAC * reference),
tolerating dips shorter than GAP bp. Genes with reference depth < MIN_REF are not called (low_expr).
Flags: hit_neighbor = extension runs to a neighbouring CDS (possible read-through from a same-strand neighbour's UTR);
hit_cap = extension reaches MAX_WIN."""
import sys, re, subprocess, statistics
from collections import defaultdict
ref, bamdir, out, samples = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:]
MAX_WIN, REF_LEN, SMOOTH, GAP = 2500, 150, 21, 30
MIN_REF, MIN_ABS, FRAC, MINQ = 10.0, 3.0, 0.15, 10

genes = defaultdict(list)
for l in open(ref):
    if l.startswith("#"): continue
    f = l.rstrip("\n").split("\t")
    if len(f) > 8 and f[2] == "mRNA":
        genes[f[0]].append([int(f[3]), int(f[4]), f[6], re.search(r"Parent=([^;]+)", f[8]).group(1)])
for c in genes: genes[c].sort()

def smooth(v, k):
    n = len(v); pre = [0.0] * (n + 1)
    for i, x in enumerate(v): pre[i + 1] = pre[i] + x
    h = k // 2
    return [(pre[min(n, i + h + 1)] - pre[max(0, i - h)]) / (min(n, i + h + 1) - max(0, i - h)) for i in range(n)]

def call(profile, ref_cov):
    """profile[i] = depth at distance i+1 from the CDS edge. returns extension length (bp)."""
    thr = max(MIN_ABS, FRAC * ref_cov)
    sm = smooth(profile, SMOOTH)
    last, gap = 0, 0
    for i, x in enumerate(sm):
        if x >= thr: last, gap = i + 1, 0
        else:
            gap += 1
            if gap >= GAP: break
    return last

def depth_array(contig, strand_label, length):
    bams = [f"{bamdir}/{s}.{strand_label}.bam" for s in samples]
    p = subprocess.Popen(["samtools", "depth", "-a", "-Q", str(MINQ), "-r", contig] + bams,
                         stdout=subprocess.PIPE, text=True)
    d = [0] * length
    for line in p.stdout:
        f = line.split("\t")
        i = int(f[1]) - 1
        if i < length: d[i] = sum(map(int, f[2:]))
    p.wait()
    return d

# contig lengths from the BAM header of the first plus bam
hdr = subprocess.run(["samtools", "view", "-H", f"{bamdir}/{samples[0]}.plus.bam"], capture_output=True, text=True).stdout
clen = {m.group(1): int(m.group(2)) for m in re.finditer(r"@SQ\tSN:(\S+)\tLN:(\d+)", hdr)}

rows = []
for c, gl in genes.items():
    L = clen[c]
    depth = {"+": depth_array(c, "plus", L), "-": depth_array(c, "minus", L)}
    sys.stderr.write(f"{c}: {len(gl)} genes\n")
    for i, (cs, ce, st, gid) in enumerate(gl):
        prev_end = gl[i - 1][1] if i > 0 else 0
        next_start = gl[i + 1][0] if i + 1 < len(gl) else L + 1
        # windows on each side of the CDS, ordered outward from the CDS edge
        left_room = cs - 1 - prev_end        # free bp left of CDS (between neighbouring CDS)
        right_room = next_start - 1 - ce
        d = depth[st]
        def outward(side):  # side 'L' or 'R' of the CDS; returns (profile, room, ref_cov)
            if side == "L":
                n = min(MAX_WIN, max(left_room, 0))
                prof = [d[cs - 2 - k] for k in range(n)]            # position cs-1-k (1-based) -> index cs-2-k
                refd = d[cs - 1: min(ce, cs - 1 + REF_LEN)]
                return prof, left_room, statistics.median(refd) if refd else 0
            n = min(MAX_WIN, max(right_room, 0))
            prof = [d[ce + k] for k in range(n)]                    # position ce+1+k -> index ce+k
            refd = d[max(cs - 1, ce - REF_LEN): ce]
            return prof, right_room, statistics.median(refd) if refd else 0
        sides = ("L", "R") if st == "+" else ("R", "L")             # 5' side first
        res = []
        for side in sides:
            prof, room, refc = outward(side)
            if refc < MIN_REF: res.append(("NA", len(prof), refc, "low_expr")); continue
            ext = call(prof, refc)
            flag = "ok"
            if ext >= len(prof) - 15 and len(prof) > 0:
                flag = "hit_cap" if room > MAX_WIN else "hit_neighbor"
            res.append((ext, len(prof), refc, flag))
        rows.append([gid, c, st, cs, ce] + [round(res[0][2], 1), round(res[1][2], 1)] + [res[0][0], res[1][0], res[0][1], res[1][1], res[0][3], res[1][3]])

with open(out, "w") as o:
    o.write("\t".join(["geneID","contig","strand","cds_start","cds_end","ref_depth5","ref_depth3","utr5_cov","utr3_cov","win5","win3","flag5","flag3"]) + "\n")
    for r in rows: o.write("\t".join(map(str, r)) + "\n")
n_ok5 = sum(1 for r in rows if r[11] == "ok" and r[7] != "NA" and r[7] > 0)
sys.stderr.write(f"genes: {len(rows)}; 5' UTR extension called (>0, flag ok): {n_ok5}\n")
