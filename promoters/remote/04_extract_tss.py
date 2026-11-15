#!/usr/bin/env python3
"""Per-gene transcript start/end (TSS/TTS) from per-sample StringTie GTFs.
usage: 04_extract_tss.py ref.gff3 out.tsv sample1.gtf sample2.gtf ...

Transcripts are assigned to genes BY COORDINATES (same contig, same strand, transcript span contains the gene's
whole CDS span), not by gene_id, so it does not matter how StringTie named the transcripts. If a transcript
contains the CDS of more than one gene (merged neighbours) it is skipped. Per sample the best-covered
transcript per gene is kept; the gene's TSS/TTS = median across samples.
Diagnostics (transcripts seen, rejected and why) are printed to stderr.
Columns: geneID contig strand cds_start cds_end tx_start tx_end n_samples tss_spread utr5_len utr3_len tss_support
(tx_start/tx_end are genomic coordinates, tx_start<=tx_end regardless of strand)."""
import sys, re, bisect, statistics
from collections import defaultdict
ref, out, gtfs = sys.argv[1], sys.argv[2], sys.argv[3:]
MIN_COV = 2.0          # StringTie per-transcript coverage

genes = {}
by_contig = defaultdict(list)
for l in open(ref):
    if l.startswith("#"): continue
    f = l.rstrip("\n").split("\t")
    if len(f) > 8 and f[2] == "mRNA":
        gid = re.search(r"Parent=([^;]+)", f[8]).group(1)
        genes[gid] = (f[0], f[6], int(f[3]), int(f[4]))      # contig, strand, cds_start, cds_end
        by_contig[f[0]].append((int(f[3]), int(f[4]), f[6], gid))
starts_idx = {}
for c, v in by_contig.items():
    v.sort(); starts_idx[c] = ([x[0] for x in v], v)

starts = {g: [] for g in genes}; ends = {g: [] for g in genes}
for gtf in gtfs:
    best = {}
    stat = defaultdict(int)
    for l in open(gtf):
        if l.startswith("#"): continue
        f = l.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "transcript": continue
        stat["transcripts"] += 1
        m = re.search(r'cov "([^"]+)"', f[8])
        cov = float(m.group(1)) if m else 0.0
        if cov < MIN_COV: stat["low_cov"] += 1; continue
        c, s, e, st = f[0], int(f[3]), int(f[4]), f[6]
        if c not in starts_idx: stat["contig_not_in_gff"] += 1; continue
        if st not in "+-": stat["no_strand"] += 1; continue
        sidx, v = starts_idx[c]
        i = bisect.bisect_left(sidx, s)
        hits = []
        while i < len(v) and v[i][0] <= e:           # genes starting inside the transcript
            if v[i][1] <= e: hits.append(v[i])        # ...and ending inside it => CDS fully contained
            i += 1
        if not hits: stat["contains_no_cds"] += 1; continue
        same = [h for h in hits if h[2] == st]
        if not same: stat["strand_mismatch"] += 1; continue
        if len(same) > 1: stat["spans_multiple_genes"] += 1; continue
        gid = same[0][3]
        stat["assigned"] += 1
        if gid not in best or cov > best[gid][0]: best[gid] = (cov, s, e)
    for gid, (cov, s, e) in best.items():
        starts[gid].append(s); ends[gid].append(e)
    sys.stderr.write(f"{gtf}: " + ", ".join(f"{k}={v}" for k, v in stat.items()) + f"; genes with a transcript: {len(best)}\n")

with open(out, "w") as o:
    o.write("\t".join(["geneID","contig","strand","cds_start","cds_end","tx_start","tx_end","n_samples","tss_spread","utr5_len","utr3_len","tss_support"]) + "\n")
    for gid, (c, st, cs, ce) in sorted(genes.items(), key=lambda x: (x[1][0], x[1][2])):
        n = len(starts[gid])
        if n:
            ts = int(statistics.median(starts[gid])); te = int(statistics.median(ends[gid]))
            five = starts[gid] if st == "+" else ends[gid]
            spread = max(five) - min(five)
            u5 = (cs - ts) if st == "+" else (te - ce)
            u3 = (te - ce) if st == "+" else (cs - ts)
        else:
            ts, te, spread, u5, u3 = cs, ce, "NA", 0, 0
        o.write("\t".join(map(str, [gid, c, st, cs, ce, ts, te, n, spread, u5, u3, "TRUE" if n else "FALSE"])) + "\n")
n_sup = sum(1 for g in genes if starts[g])
sys.stderr.write(f"genes with RNA-seq TSS support: {n_sup} of {len(genes)}\n")
