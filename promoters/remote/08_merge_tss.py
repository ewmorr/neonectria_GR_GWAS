#!/usr/bin/env python3
"""Combine StringTie TSS/TTS (04_extract_tss.py) with coverage-based UTR calls (07_utr_from_coverage.py).
usage: 08_merge_tss.py tss_table.tsv utr_coverage.tsv out.tsv
Rule per end: use StringTie when it extended the UTR; otherwise use the coverage call when it is an unambiguous
extension (flag ok, > 0); otherwise the start/stop codon. Same columns as tss_table.tsv plus utr5_source/utr3_source
(stringtie | coverage | none). Concordance of the two methods is printed for genes where both extend."""
import sys, csv, statistics
tss = {r["geneID"]: r for r in csv.DictReader(open(sys.argv[1]), delimiter="\t")}
cov = {r["geneID"]: r for r in csv.DictReader(open(sys.argv[2]), delimiter="\t")}
MIN_EXT = 20   # ignore trivial extensions
out = []; d5 = []; d3 = []; src = {"stringtie": 0, "coverage": 0, "none": 0}
for gid, t in tss.items():
    st, cs, ce = t["strand"], int(t["cds_start"]), int(t["cds_end"])
    sup = t["tss_support"] == "TRUE"
    c = cov.get(gid)
    ext = {}; source = {}
    for end, tcol, ccol, fcol in (("5", "utr5_len", "utr5_cov", "flag5"), ("3", "utr3_len", "utr3_cov", "flag3")):
        s_len = int(t[tcol]) if sup else 0
        c_len = int(c[ccol]) if c and c[ccol] != "NA" and c[fcol] == "ok" else 0
        if s_len >= MIN_EXT:
            ext[end], source[end] = s_len, "stringtie"
            if c_len >= MIN_EXT: (d5 if end == "5" else d3).append((s_len, c_len))
        elif c_len >= MIN_EXT:
            ext[end], source[end] = c_len, "coverage"
        else:
            ext[end], source[end] = 0, "none"
        src[source[end]] += 1
    left = ext["5"] if st == "+" else ext["3"]; right = ext["3"] if st == "+" else ext["5"]
    covered = sup or (c is not None and c["ref_depth5"] != "NA" and float(c["ref_depth5"]) >= 10)
    out.append([gid, t["contig"], st, cs, ce, cs - left, ce + right, t["n_samples"], t["tss_spread"], ext["5"], ext["3"],
                "TRUE" if covered else "FALSE", source["5"], source["3"]])
with open(sys.argv[3], "w") as o:
    o.write("\t".join(["geneID","contig","strand","cds_start","cds_end","tx_start","tx_end","n_samples","tss_spread","utr5_len","utr3_len","tss_support","utr5_source","utr3_source"]) + "\n")
    for r in out: o.write("\t".join(map(str, r)) + "\n")
n5 = sum(1 for r in out if r[9] > 0); n3 = sum(1 for r in out if r[10] > 0)
sys.stderr.write(f"genes {len(out)}; with 5' UTR {n5}; with 3' UTR {n3}; end sources: {src}\n")
for lab, d in (("5'", d5), ("3'", d3)):
    if len(d) > 5:
        a = [x[0] for x in d]; b = [x[1] for x in d]
        diffs = [abs(x - y) for x, y in d]
        rho = statistics.correlation(a, b) if len(d) > 2 else float("nan")
        sys.stderr.write(f"{lab} concordance, genes extended by both (n={len(d)}): median |diff| {statistics.median(diffs):.0f} bp, Pearson r {rho:.2f}\n")
