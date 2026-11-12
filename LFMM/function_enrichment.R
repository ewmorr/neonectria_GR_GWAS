## Functional enrichment of genes carrying missense / LoF-type significant SNPs (MAT, HDD4, ...)
## Terms tested: GO, PFAM, InterPro, CAZyme family, plus secreted / protease flags.
## One-sided Fisher exact test (over-representation), BH-adjusted within covariate x term type.
## Background (primary) = genes with >=1 missense/LoF-type-capable coding SNP in the dataset
## (genes without any variable coding site cannot be hits). Background "all" = all protein-coding genes.
## Run from repo root after snp_effects.R. Output: data/LFMM/function_enrichment/<cov>.enrichment.txt
library(data.table)

covariates = c("mat", "hdd4_growing")
classes = c("missense", "LoF/other_coding")
min_in_set = 2

ann = fread("data/genotype_dat/Neonectria_faginata.annotations.txt", sep = "\t", quote = "",
            select = c("GeneID", "TranscriptID", "Feature", "PFAM", "InterPro", "GO Terms", "Secreted", "Protease", "CAZyme"))
ann = ann[Feature == "mRNA"]
all_genes = ann$GeneID

# long term tables: gene, type, term, label
split_terms = function(col, type, pat = NULL){
    d = ann[, .(GeneID, v = get(col))][v != ""]
    d = d[, .(item = unlist(strsplit(v, ";", fixed = TRUE))), GeneID]
    d[, item := trimws(item)]
    d
}
go = split_terms("GO Terms", "GO"); go[, term := sub(".*(GO:[0-9]+).*", "\\1", item)]
go[, label := sub(" \\[Evidence.*", "", sub("^[^:]+: ", "", item))]; go[, type := "GO"]
pf = split_terms("PFAM", "PFAM"); pf[, `:=`(term = item, label = item, type = "PFAM")]
ip = split_terms("InterPro", "InterPro"); ip = ip[grepl("^IPR", item)]
ip[, `:=`(term = sub(" .*", "", item), label = item, type = "InterPro")]
cz = split_terms("CAZyme", "CAZyme"); cz[, `:=`(term = item, label = item, type = "CAZyme")]
flags = rbind(
    ann[Secreted != "", .(GeneID, term = "secreted(SignalP)", label = "secreted (SignalP)", type = "flag")],
    ann[Protease != "", .(GeneID, term = "protease", label = "protease (MEROPS)", type = "flag")])
terms = rbind(go[, .(GeneID, type, term, label)], pf[, .(GeneID, type, term, label)],
              ip[, .(GeneID, type, term, label)], cz[, .(GeneID, type, term, label)], flags)
terms = unique(terms)
term_lab = unique(terms[, .(type, term, label)], by = c("type", "term"))

# genes with at least one coding SNP (any class csq calls missense/synonymous/LoF) = background B
eff = fread("data/LFMM/snp_effects/all_snps.effects.txt")
bg_coding = unique(eff[effect_class %in% c("missense", "synonymous", "LoF/other_coding"), geneID])
message("background genes: coding-SNP ", length(bg_coding), "; all ", length(all_genes))

dir.create("data/LFMM/function_enrichment", showWarnings = FALSE)
run_test = function(set, bg){
    set = intersect(set, bg)
    tb = terms[GeneID %in% bg]
    n_bg = length(bg); n_set = length(set)
    cnt = tb[, .(bg_n = uniqueN(GeneID), set_n = uniqueN(GeneID[GeneID %in% set])), .(type, term)]
    cnt = cnt[set_n >= min_in_set]
    cnt[, p := mapply(function(a, b) fisher.test(matrix(c(a, n_set - a, b - a, n_bg - n_set - (b - a)), 2),
                                                 alternative = "greater")$p.value, set_n, bg_n)]
    cnt[, fdr := p.adjust(p, "BH"), by = type]
    cnt[, `:=`(n_set = n_set, n_bg = n_bg, frac_set = set_n / n_set, frac_bg = bg_n / n_bg)]
    merge(cnt, term_lab, by = c("type", "term"))[order(p)]
}

for(cov in covariates){
    s = fread(paste0("data/LFMM/annotations/", cov, ".SNPs.gene_found.nearest_neighbors.txt"))
    set = unique(unlist(strsplit(s[effect_class %in% classes & match.type == "in.gene", geneID], ";")))
    r1 = run_test(set, bg_coding); r1[, background := "coding_SNP_genes"]
    r2 = run_test(set, all_genes); r2[, background := "all_genes"]
    r = rbind(r1, r2)
    fwrite(r, paste0("data/LFMM/function_enrichment/", cov, ".enrichment.txt"), sep = "\t")
    cat("\n====", cov, ": ", length(set), "genes ====\n")
    print(r[background == "coding_SNP_genes" & p < 0.05, .(type, label = substr(label, 1, 55), set_n, bg_n, frac_set = round(frac_set, 3), frac_bg = round(frac_bg, 3), p = signif(p, 2), fdr = signif(fdr, 2))][1:min(.N, 25)])
}
