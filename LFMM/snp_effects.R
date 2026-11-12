## Parse bcftools csq output (LFMM/snp_effects.sh) into one row per SNP, then
##  1. write data/LFMM/snp_effects/all_snps.effects.txt
##  2. add effect columns to each covariate's data/LFMM/annotations/<cov>.SNPs.gene_found.nearest_neighbors.txt
##  3. compare effect classes of significant SNPs against the genome-wide background
## Run from repo root.
library(data.table)

covariates = c("topt", "t50low", "ctmax", "grpc1", "grpc2", "mat", "hdd4_growing")
out_dir = "data/LFMM/snp_effects"

# Consequence severity (most severe first); csq lists can include several (e.g. splice_region + missense)
sev = c("start_lost", "stop_gained", "stop_lost", "splice_acceptor", "splice_donor", "frameshift",
        "inframe_altering", "missense", "splice_region", "start_retained", "stop_retained",
        "synonymous", "intron", "non_coding", "inframe_insertion", "inframe_deletion")
# collapsed class used for comparisons
cls = function(x) fcase(
    x %in% c("missense"), "missense",
    x %in% c("synonymous", "stop_retained", "start_retained"), "synonymous",
    x %in% c("stop_gained", "start_lost", "stop_lost", "splice_acceptor", "splice_donor", "frameshift", "inframe_altering"), "LoF/other_coding",
    x %in% c("intron", "splice_region"), "intron/splice_region",
    x == "intergenic", "intergenic",
    default = "other")

x = fread(file.path(out_dir, "all_snps.csq.txt"), header = FALSE, sep = "\t",
          col.names = c("scaffold", "position", "ref", "alt", "bcsq"))

parse_one = function(b){
    if(b == "." || is.na(b)) return(list(effect = "intergenic", geneID = NA_character_, aa_change = NA_character_, n_conseq = 0L))
    e = strsplit(b, ",", fixed = TRUE)[[1]]
    e = sub("^[@*]", "", e)                      # drop haplotype markers if any
    f = strsplit(e, "|", fixed = TRUE)
    type = sapply(f, `[`, 1)
    # a single csq term can be "a&b"; rank by the most severe component
    rk = sapply(strsplit(type, "&", fixed = TRUE), function(t) min(match(t, sev, nomatch = length(sev) + 1)))
    k = which.min(rk)
    top = strsplit(type[k], "&", fixed = TRUE)[[1]]
    top = top[which.min(match(top, sev, nomatch = length(sev) + 1))]
    list(effect = top, geneID = f[[k]][2], aa_change = if(length(f[[k]]) >= 6) f[[k]][6] else NA_character_, n_conseq = length(e))
}
p = rbindlist(lapply(x$bcsq, parse_one))
x = cbind(x[, .(scaffold, position, ref, alt)], p)
x[, effect_class := cls(effect)]
fwrite(x, file.path(out_dir, "all_snps.effects.txt"), sep = "\t")
cat("All SNPs by class:\n"); print(x[, .N, effect_class][order(-N)])

# --- per covariate ---
summ = list()
for(cov in covariates){
    f = paste0("data/LFMM/annotations/", cov, ".SNPs.gene_found.nearest_neighbors.txt")
    s = fread(f)
    s[, c("effect", "effect_class", "aa_change", "ref", "alt") := NULL][]  # allow re-running
    s = merge(s, x[, .(scaffold, position, ref, alt, effect, effect_class, aa_change)],
              by = c("scaffold", "position"), all.x = TRUE, sort = FALSE)
    stopifnot(!anyNA(s$effect_class) || nrow(s) == 0)
    fwrite(s, f, sep = "\t")
    if(nrow(s) > 0) summ[[cov]] = s[, .(cov = cov, N = .N), effect_class]
}

# --- sig vs background enrichment (Fisher) ---
bg = x[, .N, effect_class]
res = rbindlist(lapply(names(summ), function(cov){
    s = summ[[cov]]
    rbindlist(lapply(bg$effect_class, function(k){
        a = s[effect_class == k, sum(N)]; nsig = sum(s$N)
        b = bg[effect_class == k, N]; nbg = sum(bg$N)
        data.table(cov, effect_class = k, n_sig = a, n_sig_total = nsig, frac_sig = a / nsig,
                   frac_all = b / nbg,
                   fisher_p = fisher.test(matrix(c(a, nsig - a, b - a, nbg - nsig - (b - a)), 2))$p.value)
    }))
}))
fwrite(res, file.path(out_dir, "sig_vs_background_effect_classes.txt"), sep = "\t")
print(res[order(cov, effect_class)], digits = 3)
