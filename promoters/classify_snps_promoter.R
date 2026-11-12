## Classify SNPs relative to gene structure: coding / 5'UTR / 3'UTR / promoter (core, proximal) / 3' flank / distal intergenic.
## Uses transcript start/end per gene (promoters/remote/04_extract_tss.py output) when present;
## otherwise falls back to CDS coordinates (the start codon stands in for the TSS: promoter distances are then
## measured from the ATG and UTR classes are empty). Run from the repo root.
##   Rscript promoters/classify_snps_promoter.R [tss_table.tsv]
## Output: data/promoters/all_snps.region.txt, per-covariate region columns appended to
##         data/LFMM/annotations/<cov>.SNPs.gene_found.nearest_neighbors.txt, and sig-vs-background tests.
library(data.table)

args = commandArgs(trailingOnly = TRUE)
tss_file = if(length(args) >= 1) args[1] else "data/promoters/tss_table.tsv"
CORE = 250       # bp upstream of TSS: "core" promoter window
PROX = 1000      # bp upstream of TSS: proximal promoter window (candidate regulatory region)
FLANK3 = 500     # bp downstream of the transcript end: 3' flank (terminator region)
covariates = c("topt", "t50low", "ctmax", "grpc1", "grpc2", "mat", "hdd4_growing")
dir.create("data/promoters", showWarnings = FALSE)

# ---------- gene table ----------
gff = fread("data/genotype_dat/Neonectria_faginata.gff3", sep = "\t", header = FALSE, skip = "tig", quote = "",
            col.names = c("contig", "src", "type", "start", "stop", "score", "strand", "phase", "attr"))
g = gff[type == "mRNA", .(contig, strand, cds_start = start, cds_end = stop,
                          geneID = sub(".*Parent=([^;]+).*", "\\1", attr))]
if(file.exists(tss_file)){
    t = fread(tss_file)[, .(geneID, tx_start, tx_end, tss_support)]
    g = merge(g, t, by = "geneID", all.x = TRUE)
    g[is.na(tss_support), tss_support := FALSE]
    message("TSS table: ", sum(g$tss_support), " of ", nrow(g), " genes with RNA-seq transcript support")
}else{
    message("No TSS table found (", tss_file, "): using CDS coordinates as transcript extents")
    g[, `:=`(tx_start = cds_start, tx_end = cds_end, tss_support = FALSE)]
}
g[is.na(tx_start), tx_start := cds_start]; g[is.na(tx_end), tx_end := cds_end]
setorder(g, contig, cds_start)
# StringTie can merge neighbouring genes in compact genomes: do not let UTRs run into a neighbour's CDS
g[, prev_cds_end := shift(cds_end, 1, fill = 0L), by = contig]
g[, next_cds_start := shift(cds_start, 1, type = "lead", fill = .Machine$integer.max), by = contig]
g[, `:=`(clipped = tx_start <= prev_cds_end | tx_end >= next_cds_start)]
g[, `:=`(tx_start = pmax(tx_start, prev_cds_end + 1L), tx_end = pmin(tx_end, next_cds_start - 1L))]
g[, `:=`(tx_start = pmin(tx_start, cds_start), tx_end = pmax(tx_end, cds_end))]
message("transcript extents clipped at a neighbouring CDS: ", sum(g$clipped), " genes")

# ---------- SNPs ----------
snp = fread("data/LFMM/snp_effects/all_snps.effects.txt")[, .(scaffold, position, effect_class)]
setnames(snp, c("scaffold", "position"), c("contig", "pos"))
snp[, id := .I]
snp[, c("region_class", "promoter_gene", "promoter_dist", "promoter_tss_support", "flank3_gene", "flank3_dist", "bidirectional") :=
        .(NA_character_, NA_character_, NA_integer_, NA, NA_character_, NA_integer_, FALSE)]

# SNPs inside a gene's transcript extent (CDS vs UTR)
ov = g[snp, on = .(contig, tx_start <= pos, tx_end >= pos), nomatch = NULL,
       .(id = i.id, geneID, strand = x.strand, pos = i.pos, cds_start = x.cds_start, cds_end = x.cds_end, tss_support = x.tss_support)]
ov[, cls := fifelse(pos >= cds_start & pos <= cds_end, "coding",
             fifelse((strand == "+" & pos < cds_start) | (strand == "-" & pos > cds_end), "5UTR", "3UTR"))]
ov = ov[order(id, cls != "coding")][!duplicated(id)]   # if overlapping extents, coding wins
snp[ov$id, `:=`(region_class = ov$cls, promoter_gene = ov$geneID)]

# intergenic SNPs: nearest transcript end to the left (L) and nearest transcript start to the right (R)
ig = snp[is.na(region_class), .(id, contig, pos)]
gl = g[, .(contig, pos = tx_end, L_gene = geneID, L_strand = strand, L_end = tx_end, L_sup = tss_support)]
gr = g[, .(contig, pos = tx_start, R_gene = geneID, R_strand = strand, R_start = tx_start, R_sup = tss_support)]
setkey(gl, contig, pos); setkey(gr, contig, pos)
ig = gl[ig, on = .(contig, pos), roll = Inf]
ig = gr[ig, on = .(contig, pos), roll = -Inf]
ig[, `:=`(dL = pos - L_end, dR = R_start - pos)]
# L on "-" strand: SNP is upstream of L (L's TSS is at its right end). R on "+" strand: SNP is upstream of R.
ig[, `:=`(L_up = !is.na(L_gene) & L_strand == "-", R_up = !is.na(R_gene) & R_strand == "+")]
ig[, `:=`(L_prom = L_up & dL <= PROX, R_prom = R_up & dR <= PROX)]
ig[, `:=`(pick_L = L_prom & (!R_prom | dL <= dR))]
ig[, `:=`(pgene = fifelse(pick_L, L_gene, fifelse(R_prom, R_gene, NA_character_)),
          pdist = fifelse(pick_L, as.integer(dL), fifelse(R_prom, as.integer(dR), NA_integer_)),
          psup  = fifelse(pick_L, L_sup, fifelse(R_prom, R_sup, NA)),
          bidir = L_prom & R_prom)]
# 3' flank: L on "+" (SNP downstream of its end) or R on "-"
ig[, `:=`(L_t = !is.na(L_gene) & L_strand == "+" & dL <= FLANK3, R_t = !is.na(R_gene) & R_strand == "-" & dR <= FLANK3)]
ig[, `:=`(tgene = fifelse(L_t & (!R_t | dL <= dR), L_gene, fifelse(R_t, R_gene, NA_character_)),
          tdist = fifelse(L_t & (!R_t | dL <= dR), as.integer(dL), fifelse(R_t, as.integer(dR), NA_integer_)))]
ig[, region_class := fifelse(!is.na(pgene) & pdist <= CORE, "promoter_core",
                     fifelse(!is.na(pgene), "promoter_proximal",
                     fifelse(!is.na(tgene), "3prime_flank", "distal_intergenic")))]
snp[ig$id, `:=`(region_class = ig$region_class, promoter_gene = ig$pgene, promoter_dist = ig$pdist,
                promoter_tss_support = ig$psup, flank3_gene = ig$tgene, flank3_dist = ig$tdist, bidirectional = ig$bidir)]
snp[region_class %in% c("coding", "5UTR", "3UTR"), promoter_dist := NA_integer_]

out = snp[, .(scaffold = contig, position = pos, region_class, promoter_gene, promoter_dist, promoter_tss_support,
              flank3_gene, flank3_dist, bidirectional)]
fwrite(out, "data/promoters/all_snps.region.txt", sep = "\t")
cat("\nAll SNPs by region class:\n"); print(out[, .(N = .N, frac = round(.N / nrow(out), 4)), region_class][order(-N)])

# ---------- significant SNPs vs background ----------
bg = out[, .N, region_class]
res = list()
for(cov in covariates){
    f = paste0("data/LFMM/annotations/", cov, ".SNPs.gene_found.nearest_neighbors.txt")
    s = fread(f)
    s[, intersect(names(s), c("region_class", "promoter_gene", "promoter_dist", "promoter_tss_support", "flank3_gene", "flank3_dist", "bidirectional")) := NULL][]
    s = merge(s, out, by = c("scaffold", "position"), all.x = TRUE, sort = FALSE)
    fwrite(s, f, sep = "\t")
    if(nrow(s) == 0) next
    for(k in bg$region_class){
        a = sum(s$region_class == k); nsig = nrow(s); b = bg[region_class == k, N]; nbg = sum(bg$N)
        res[[length(res) + 1]] = data.table(cov, region_class = k, n_sig = a, n_sig_total = nsig, frac_sig = a / nsig,
            frac_all = b / nbg, fisher_p = fisher.test(matrix(c(a, nsig - a, b - a, nbg - nsig - (b - a)), 2))$p.value)
    }
}
res = rbindlist(res)
fwrite(res, "data/promoters/sig_vs_background_region_classes.txt", sep = "\t")
print(res[order(cov, region_class)], digits = 3)
