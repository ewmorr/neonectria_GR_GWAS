#!/usr/bin/env python3
"""Per-gene transcript start/end (TSS/TTS) from per-sample StringTie GTFs.
usage: 04_extract_tss.py ref.gff3 out.tsv sample1.gtf sample2.gtf ...
For every reference gene: in each sample take the best-covered transcript carrying that gene_id that spans the
CDS start codon (so it is the 5'-complete isoform); the gene's TSS/TTS = median across samples.
Columns: geneID contig strand cds_start cds_end tx_start tx_end n_samples tss_spread utr5_len utr3_len tss_support
(tx_start/tx_end are in genomic coordinates, tx_start<=tx_end regardless of strand)."""
import sys, re, statistics
ref, out, gtfs = sys.argv[1], sys.argv[2], sys.argv[3:]
MIN_COV = 2.0          # StringTie per-transcript coverage
genes = {}
for l in open(ref):
    if l.startswith("#"): continue
    f = l.rstrip("\n").split("\t")
    if len(f) > 8 and f[2] == "mRNA":
        gid = re.search(r"Parent=([^;]+)", f[8]).group(1)
        genes[gid] = (f[0], f[6], int(f[3]), int(f[4]))      # contig, strand, cds_start(=gene start), cds_end
starts = {g: [] for g in genes}; ends = {g: [] for g in genes}
for gtf in gtfs:
    best = {}
    for l in open(gtf):
        if l.startswith("#"): continue
        f = l.rstrip("\n").split("\t")
        if f[2] != "transcript": continue
        gid = re.search(r'gene_id "([^"]+)"', f[8]).group(1)
        if gid not in genes: continue
        cov = float(re.search(r'cov "([^"]+)"', f[8]).group(1))
        s, e = int(f[3]), int(f[4])
        c, st, cs, ce = genes[gid]
        if cov < MIN_COV or not (s <= cs and e >= ce): continue   # must contain the whole CDS span
        if gid not in best or cov > best[gid][0]: best[gid] = (cov, s, e)
    for gid, (cov, s, e) in best.items():
        starts[gid].append(s); ends[gid].append(e)
with open(out, "w") as o:
    o.write("\t".join(["geneID","contig","strand","cds_start","cds_end","tx_start","tx_end","n_samples","tss_spread","utr5_len","utr3_len","tss_support"]) + "\n")
    for gid, (c, st, cs, ce) in sorted(genes.items(), key=lambda x: (x[1][0], x[1][2])):
        n = len(starts[gid])
        if n:
            ts = int(statistics.median(starts[gid])); te = int(statistics.median(ends[gid]))
            # spread of the 5'-most end across samples (bp): relevant end depends on strand
            five = starts[gid] if st == "+" else ends[gid]
            spread = max(five) - min(five)
            u5 = (cs - ts) if st == "+" else (te - ce)
            u3 = (te - ce) if st == "+" else (cs - ts)
        else:
            ts, te, spread, u5, u3 = cs, ce, "NA", 0, 0
        o.write("\t".join(map(str, [gid, c, st, cs, ce, ts, te, n, spread, u5, u3, "TRUE" if n else "FALSE"])) + "\n")
