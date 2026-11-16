## Shared loader for the eggNOG-mapper annotations (renamed to Q2T73_ transcript IDs) with readable KEGG names.
## source("LFMM/annotation_helpers.R"); em = load_emapper()   # one row per gene (GeneID)
library(data.table)

# KEGG global/overview maps carry no functional information for enrichment
KEGG_GLOBAL = c("01100", "01110", "01120", "01130", "01200", "01210", "01212", "01230", "01232", "01240", "01250")

load_emapper = function(
        emapper = "data/genotype_dat/Nf_.emapper.annotations.renamed",
        ann = "data/genotype_dat/Neonectria_faginata.annotations.txt",
        kegg_dir = "data/genotype_dat/kegg"){
    e = fread(cmd = paste0("grep -v '^##' ", emapper), sep = "\t", quote = "", na.strings = "-")
    setnames(e, "#query", "TranscriptID")
    a = fread(ann, sep = "\t", quote = "", select = c("GeneID", "TranscriptID", "Feature"))[Feature == "mRNA", .(GeneID, TranscriptID)]
    e = merge(a, e, by = "TranscriptID")
    pw = fread(file.path(kegg_dir, "kegg_pathway_names.txt"), header = FALSE, sep = "\t", col.names = c("id", "name"))
    pw[, num := sub("^[a-z]+", "", id)]; pw = unique(pw, by = "num")
    md = fread(file.path(kegg_dir, "kegg_module_names.txt"), header = FALSE, sep = "\t", col.names = c("id", "name"))
    pwname = setNames(pw$name, pw$num); mdname = setNames(md$name, md$id)
    pw_ids = function(x) { v = unique(sub("^ko", "", grep("^ko", strsplit(x, ",", fixed = TRUE)[[1]], value = TRUE))); setdiff(v, KEGG_GLOBAL) }
    e[, kegg_pathway_ids := vapply(KEGG_Pathway, function(x) if(is.na(x)) NA_character_ else { v = pw_ids(x); if(length(v)) paste(v, collapse = ",") else NA_character_ }, "")]
    e[, kegg_pathways := vapply(kegg_pathway_ids, function(x) if(is.na(x)) NA_character_ else
        paste(paste0(strsplit(x, ",")[[1]], " ", pwname[strsplit(x, ",")[[1]]]), collapse = "; "), "")]
    e[, kegg_modules := vapply(KEGG_Module, function(x) if(is.na(x)) NA_character_ else
        paste(paste0(strsplit(x, ",")[[1]], " ", mdname[strsplit(x, ",")[[1]]]), collapse = "; "), "")]
    e[, .(geneID = GeneID, eggnog_name = Preferred_name, eggnog_description = Description, COG_category,
          EC_eggnog = EC, KEGG_ko = sub("ko:", "", KEGG_ko), kegg_pathway_ids, kegg_pathways, KEGG_Module, kegg_modules,
          GOs_eggnog = GOs, BRITE)]
}
