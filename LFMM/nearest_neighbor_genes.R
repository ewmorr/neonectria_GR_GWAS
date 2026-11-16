## Map SNPs significantly associated (LFMM) with each covariate to genes:
## genes containing the SNP ("in.gene"), otherwise the nearest protein-coding gene(s) ("nearest").
## Run from the repo root. Replaces the per-covariate nearest_neighbor_genes.*.r scripts.
##
## Distance sign: genes downstream of the SNP (higher coordinate) are positive, upstream negative
## (gene start - SNP pos, or gene stop - SNP pos). Ties / overlapping genes are all kept; geneID
## holds them separated by ";" and each gene gets its own SNP count.

library(dplyr)
library(data.table)
library(ggplot2)
source("LFMM/annotation_helpers.R")

covariates = c("topt", "t50low", "ctmax", "grpc1", "grpc2", "mat", "hdd4_growing")

############################
# Gene models (protein-coding only = mRNA features)
############################
gff = fread("data/genotype_dat/Neonectria_faginata.gff3", sep = "\t", header = FALSE,
            skip = "tig", quote = "", col.names = c("contig", "source", "type", "start", "stop",
            "score", "strand", "phase", "attr"))
genes = gff[type == "mRNA"]
genes[, geneID := sub(".*Parent=([^;]+).*", "\\1", attr)]
genes[, transcriptID := sub(".*ID=([^;]+).*", "\\1", attr)]
genes = genes[, .(contig, geneID, transcriptID, start, stop, strand)]
stopifnot(!anyDuplicated(genes$geneID))

# Functional annotation (keyed by TranscriptID); drop the big sequence columns
ann_cols = c("TranscriptID", "Name", "Product", "EC_number", "PFAM", "InterPro", "EggNog",
             "COG", "GO Terms", "Secreted", "Membrane", "Protease", "CAZyme", "antiSMASH")
ann = fread("data/genotype_dat/Neonectria_faginata.annotations.txt", sep = "\t",
            select = ann_cols, quote = "")
setnames(ann, "TranscriptID", "transcriptID")
genes = merge(genes, ann, by = "transcriptID", all.x = TRUE, sort = FALSE)
em = load_emapper()   # eggNOG-mapper: description, KEGG ko/pathways/modules, COG, EC
genes = merge(genes, em, by = "geneID", all.x = TRUE, sort = FALSE)
setkey(genes, contig, start)

############################
# Per-SNP matching
############################
# Returns list(type, idx, dist) for one SNP; idx = row(s) of `g` (genes on the SNP's contig)
match_snp = function(pos, g){
    inside = which(g$start <= pos & g$stop >= pos)
    if(length(inside) > 0) return(list(type = "in.gene", idx = inside, dist = 0))
    if(nrow(g) == 0) return(list(type = "no.match", idx = integer(0), dist = NA_real_))
    d = ifelse(g$start > pos, g$start - pos, g$stop - pos)
    nearest = which(abs(d) == min(abs(d)))
    list(type = "nearest", idx = nearest, dist = d[nearest[1]])
}

dir.create("data/LFMM/annotations", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/LFMM", showWarnings = FALSE, recursive = TRUE)

for(cov in covariates){
    message("== ", cov)
    res = read.table(paste0("data/LFMM/", cov, ".lfmm.txt"), header = TRUE)
    # manually adjusted FDR significance (as in the previous scripts)
    sig = res %>% filter(FDR.sig.man == "sig") %>% select(scaffold, position)
    sig$scaffold = sub("_pilon", "", sig$scaffold)
    message(nrow(sig), " significant SNPs")

    out_prefix = paste0("data/LFMM/annotations/", cov)
    gene_counts = copy(genes)
    gene_counts[, sig_SNP_count := 0L]
    sig$match.type = character(nrow(sig))
    sig$geneID = character(nrow(sig))
    sig$distance = numeric(nrow(sig))

    for(i in seq_len(nrow(sig))){
        g = genes[contig == sig$scaffold[i]]
        m = match_snp(sig$position[i], g)
        sig$match.type[i] = m$type
        sig$distance[i] = m$dist
        if(length(m$idx) > 0){
            ids = g$geneID[m$idx]
            sig$geneID[i] = paste(ids, collapse = ";")
            gene_counts[geneID %in% ids, sig_SNP_count := sig_SNP_count + 1L]
        }else{
            sig$geneID[i] = "no.match"
        }
    }

    hits = gene_counts[sig_SNP_count > 0][order(-sig_SNP_count)]
    if(nrow(sig) > 0){
        message("match types: ", paste(names(table(sig$match.type)), table(sig$match.type), collapse = ", "))
        message("genes with SNPs: ", nrow(hits), "; genes with >1 SNP: ", sum(hits$sig_SNP_count > 1))
        nn = abs(sig$distance[sig$match.type == "nearest"])
        if(length(nn) > 0) message("nearest-gene |distance| median ", median(nn), ", range ", min(nn), "-", max(nn))
    }

    # tables (written even if empty so downstream files reflect the current annotation)
    write.table(hits, paste0(out_prefix, ".gene_SNP_hits.txt"), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
    write.table(sig, paste0(out_prefix, ".SNPs.gene_found.nearest_neighbors.txt"), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
    write.table(hits[, .(geneID)], paste0(out_prefix, ".geneIDs.nearest_neighbors.txt"), col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)

    # plots
    if(nrow(sig) > 0){
        p1 = ggplot(sig, aes(x = distance)) +
            geom_histogram(breaks = seq(-15000, 15000, 250)) +
            theme_bw() +
            labs(x = "Gene distance from SNP (bp)", y = "Number SNPs (bin width 250 bp)")
        ggsave(paste0("figures/LFMM/SNP_distance.", cov, ".pdf"), p1, width = 6, height = 4)

        p2 = ggplot(hits, aes(x = sig_SNP_count)) +
            geom_bar() +
            labs(x = "SNPs per gene", y = "No. genes") +
            theme_bw()
        ggsave(paste0("figures/LFMM/SNP_per_gene.", cov, ".pdf"), p2, width = 6, height = 4)
    }
}
