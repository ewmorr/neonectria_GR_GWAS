# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

Collection of standalone R scripts (no build/lint/test system) analyzing growth-rate (GR) thermal performance of *Neonectria faginata* (Nf) and *N. ditissima* (Nd) isolates, and testing for genetic signal (GWAS/LFMM, selective sweeps) in the resulting parameters. It is an RStudio project (`neonectria_GR_GWAS.Rproj`).

- **Run scripts from the repo root** — scripts use root-relative paths (`data/...`, `source("library/ggplot_theme.txt")`, `source("growth_curves/make_site_metadata.r")`).
- `data/` and `figures/` are gitignored; data is not in the repo. Don't expect to run scripts end-to-end without it.
- Files are mixed `.R`/`.r` extensions; the README refers to some with lowercase names.

## Pipeline (in order)

1. `growth_curves/format_tables_calc_gr.R` — reads raw plate-reader CSVs from `data/processed_data/` (filenames encode set/rep/temperature, e.g. `set1_rep2_25C.csv`, plus repeat runs in `redos/`) and computes growth rates.
2. `growth_curves/rTPC/` — fits thermal performance curves with the `rTPC` package (`rTPC.multiple_curves.R` compares model fits; outputs like `data/rTPC/all_fits.params.csv`, `sharpeschoolhigh_half_max_temps.csv`, `d_ic.bestTwo_fits.fit_stats.csv`). Parameters used downstream: rmax/GRmax, topt, ctmax (sharpeschoolhigh and boatman), t50low/t50high, breadth, skewness. Best-model selection is done via information-criteria weights (`sharpeschoolhigh` and `boatman` are the models carried forward).
   - `spp_comparisons.R` (Nf vs Nd), `local_adaptation_lms.R` (parameters vs. climate at collection site; species interaction matters), `parameter_correlations.R`, `parameter_distributions.R`.
3. `LFMM/` — latent-factor mixed models (`lfmm` package) of Topt, t50low, ctmax against genotype data (`LFMM.topt-t50low-ctmax.r`); `nearest_neighbor_genes.R` loops over all covariates and maps significant SNPs to genes (in-gene, else nearest protein-coding gene) using `data/genotype_dat/Neonectria_faginata.gff3` + `.annotations.txt`, writing to `data/LFMM/annotations/`; `snp_effects.sh` (needs `conda activate bcftools`; converts the gff to Ensembl-style for `bcftools csq -l`) + `snp_effects.R` classify every SNP as synonymous/missense/etc. and append effect columns to those tables; `promoters/` holds the SLURM RNA-seq/StringTie TSS pipeline (`remote/`, see its README) and `classify_snps_promoter.R` (SNP region classes: promoter/UTR/flank); `function_enrichment.R` tests GO/PFAM/InterPro/CAZyme enrichment of missense/LoF gene sets; `cluster_response_vars.R` clusters the response variables. `sweepfinder2_*` are SweepFinder2 selective-sweep scripts (formatting shell script + SLURM array job run on an HPC cluster, one task per contig `tig*.sf2`, not locally).
4. `growth_curves/make_site_metadata.r` — shared site/climate metadata, `source()`d by several analysis scripts.

Other directories: `sample_randomization/` (isolate selection for growth experiments), `site_maps_and_climate/` (maps/climate panels), `R_scripts/` (misc), `growth_curves/old_param_ests/` (superseded nls-based estimation; ignore unless needed), `library/ggplot_theme.txt` (shared ggplot theme, sourced not run).

## Conventions

- Tidyverse style (`dplyr`, `ggplot2`, `lubridate`); `.Rhistory`/`.RData` are tracked/present at the root — avoid committing changes to them.
- Research notes and key findings (species differences, climate correlations) are kept in `README.md`.
