## Functional enrichment of gene sets derived from significant LFMM SNPs.
## Terms: funannotate GO / PFAM / InterPro / CAZyme / secreted / protease, plus eggNOG-mapper KEGG pathway,
## KEGG module, KEGG ortholog, COG category, EC and eggNOG GO (data/genotype_dat/Nf_.emapper.annotations.renamed).
## One-sided Fisher exact test (over-representation), BH-adjusted within set x term type; terms need >= min_in_set genes.
## Gene sets:
##   <cov>.codingLoF : genes with an in-gene missense/LoF-type significant SNP   (bg = genes with any coding SNP; also all genes)
##   <cov>.allhits   : every gene assigned to a significant SNP (in-gene OR nearest)  (bg = all genes; nearest-gene
##                     assignment is weak for intergenic SNPs, so treat as exploratory)
##   growth.allhits  : pooled topt + ctmax + grpc1 + grpc2 (each too small alone)
## Run from repo root after nearest_neighbor_genes.R and snp_effects.R.  Output: data/LFMM/function_enrichment/<set>.enrichment.txt
library(data.table)
source("LFMM/annotation_helpers.R")

classes = c("missense", "LoF/other_coding")
min_in_set = 2
growth = c("topt", "ctmax", "grpc1", "grpc2")

ann = fread("data/genotype_dat/Neonectria_faginata.annotations.txt", sep = "\t", quote = "",
            select = c("GeneID", "Feature", "PFAM", "InterPro", "GO Terms", "Secreted", "Protease", "CAZyme"))
ann = ann[Feature == "mRNA"]
all_genes = ann$GeneID
em = load_emapper()

split_terms = function(d, col){
    d = d[, .(GeneID = geneID, v = get(col))][!is.na(v) & v != ""]
    d[, .(item = trimws(unlist(strsplit(v, "[;,]")))), GeneID][item != ""]
}
ann[, geneID := GeneID]
go = split_terms(ann, "GO Terms"); go[, `:=`(term = sub(".*(GO:[0-9]+).*", "\\1", item),
      label = sub(" \\[Evidence.*", "", sub("^[^:]+: ", "", item)), type = "GO")]
pf = split_terms(ann, "PFAM");     pf[, `:=`(term = item, label = item, type = "PFAM")]
ipd = ann[InterPro != "", .(item = trimws(unlist(strsplit(InterPro, ";")))), by = GeneID][grepl("^IPR", item)]
ipd[, `:=`(term = sub(" .*", "", item), label = item, type = "InterPro")]
cz = split_terms(ann, "CAZyme");   cz[, `:=`(term = item, label = item, type = "CAZyme")]
flags = rbind(ann[Secreted != "", .(GeneID, term = "secreted", label = "secreted (SignalP)", type = "flag")],
              ann[Protease != "", .(GeneID, term = "protease", label = "protease (MEROPS)", type = "flag")])

# eggNOG-mapper terms
pwn = em[!is.na(kegg_pathway_ids), .(item = unlist(strsplit(kegg_pathways, "; "))), by = geneID]
pwn[, `:=`(term = sub(" .*", "", item), label = item, type = "KEGG_pathway")]
mdn = em[!is.na(kegg_modules), .(item = unlist(strsplit(kegg_modules, "; "))), by = geneID]
mdn[, `:=`(term = sub(" .*", "", item), label = item, type = "KEGG_module")]
kon = em[!is.na(KEGG_ko), .(item = unlist(strsplit(KEGG_ko, ","))), by = geneID]
kon[, `:=`(term = item, label = item, type = "KEGG_ko")]
cog_names = c(J = "Translation", A = "RNA processing", K = "Transcription", L = "Replication/repair", B = "Chromatin",
  D = "Cell cycle", Y = "Nuclear structure", V = "Defense", T = "Signal transduction", M = "Cell wall/membrane", N = "Motility",
  Z = "Cytoskeleton", W = "Extracellular", U = "Trafficking/secretion", O = "Chaperones/PTM", C = "Energy", G = "Carbohydrate metabolism",
  E = "Amino acid metabolism", F = "Nucleotide metabolism", H = "Coenzyme metabolism", I = "Lipid metabolism", P = "Inorganic ion transport",
  Q = "Secondary metabolites", R = "General function only", S = "Function unknown")
cgn = em[!is.na(COG_category), .(item = strsplit(COG_category, "")[[1]]), by = geneID][item %in% names(cog_names)]
cgn[, `:=`(term = item, label = paste0(item, " ", cog_names[item]), type = "COG")]
ecn = em[!is.na(EC_eggnog), .(item = unlist(strsplit(EC_eggnog, ","))), by = geneID]
ecn[, `:=`(term = item, label = item, type = "EC")]
gen = em[!is.na(GOs_eggnog), .(item = unlist(strsplit(GOs_eggnog, ","))), by = geneID]
gen[, `:=`(term = item, label = item, type = "GO_eggnog")]

cols = c("geneID", "type", "term", "label")
setnames(go, "GeneID", "geneID"); setnames(pf, "GeneID", "geneID"); setnames(cz, "GeneID", "geneID")
setnames(ipd, "GeneID", "geneID"); setnames(flags, "GeneID", "geneID")
terms = unique(rbind(go[, ..cols], pf[, ..cols], ipd[, ..cols], cz[, ..cols], flags[, ..cols],
                     pwn[, ..cols], mdn[, ..cols], kon[, ..cols], cgn[, ..cols], ecn[, ..cols], gen[, ..cols]))
term_lab = unique(terms[, .(type, term, label)], by = c("type", "term"))
# GO_eggnog labels: use GO IDs only (names not shipped with emapper); funannotate GO carries the names when overlapping
eff = fread("data/LFMM/snp_effects/all_snps.effects.txt")
bg_coding = unique(eff[effect_class %in% c("missense", "synonymous", "LoF/other_coding"), geneID])
message("annotated genes with eggNOG: ", nrow(em), " of ", length(all_genes), "; coding-SNP genes: ", length(bg_coding))

run_test = function(set, bg){
    set = intersect(set, bg)
    tb = terms[geneID %in% bg]
    n_bg = length(bg); n_set = length(set)
    cnt = tb[, .(bg_n = uniqueN(geneID), set_n = uniqueN(geneID[geneID %in% set])), .(type, term)][set_n >= min_in_set]
    if(nrow(cnt) == 0) return(NULL)
    cnt[, p := mapply(function(a, b) fisher.test(matrix(c(a, n_set - a, b - a, n_bg - n_set - (b - a)), 2), alternative = "greater")$p.value, set_n, bg_n)]
    cnt[, fdr := p.adjust(p, "BH"), by = type]
    cnt[, `:=`(n_set = n_set, n_bg = n_bg, frac_set = set_n / n_set, frac_bg = bg_n / n_bg)]
    merge(cnt, term_lab, by = c("type", "term"))[order(p)]
}
read_sig = function(cov) fread(paste0("data/LFMM/annotations/", cov, ".SNPs.gene_found.nearest_neighbors.txt"))
genes_of = function(x) unique(unlist(strsplit(x[geneID != "no.match" & geneID != "", geneID], ";")))

sets = list()
for(cov in c("mat", "hdd4_growing")){
    s = read_sig(cov)
    sets[[paste0(cov, ".codingLoF")]] = list(genes = genes_of(s[effect_class %in% classes & match.type == "in.gene"]), bg = list(coding_SNP_genes = bg_coding, all_genes = all_genes))
    sets[[paste0(cov, ".allhits")]] = list(genes = genes_of(s), bg = list(all_genes = all_genes))
}
gs = rbindlist(lapply(growth, read_sig))
sets[["growth.allhits"]] = list(genes = genes_of(gs), bg = list(all_genes = all_genes))
sets[["growth.codingLoF"]] = list(genes = genes_of(gs[effect_class %in% classes & match.type == "in.gene"]), bg = list(coding_SNP_genes = bg_coding))

dir.create("data/LFMM/function_enrichment", showWarnings = FALSE)
unlink(list.files("data/LFMM/function_enrichment", pattern = "^(mat|hdd4_growing)\\.enrichment\\.txt$", full.names = TRUE))
for(nm in names(sets)){
    r = rbindlist(lapply(names(sets[[nm]]$bg), function(b){ x = run_test(sets[[nm]]$genes, sets[[nm]]$bg[[b]]); if(!is.null(x)) x[, background := b]; x }), fill = TRUE)
    fwrite(r, paste0("data/LFMM/function_enrichment/", nm, ".enrichment.txt"), sep = "\t")
    cat("\n====", nm, ":", length(sets[[nm]]$genes), "genes ====\n")
    if(nrow(r) > 0){
        r1 = r[background == names(sets[[nm]]$bg)[1] & fdr < 0.1]
        if(nrow(r1) > 0) print(r1[order(fdr), .(type, label = substr(label, 1, 60), set_n, bg_n, p = signif(p, 2), fdr = signif(fdr, 2))][1:min(.N, 20)])
        else cat("nothing at FDR < 0.1; best p:", signif(min(r$p), 2), "\n")
    }
}
