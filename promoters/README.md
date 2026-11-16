# Promoter / TSS pipeline

Goal: classify SNPs (esp. intergenic ones) as promoter-associated using RNA-seq-supported transcript starts.

## On the cluster (SLURM) — `promoters/remote/`
1. `bash 00_env_setup.sh` once (conda env `promoter`: stringtie, samtools, hisat2, meme). Copy `data/genotype_dat/Neonectria_faginata.{scaffolds.fa,gff3}` to `$WORKDIR/ref/`; copy `samples.template.tsv` to `$WORKDIR/samples.tsv` and fill in the 8 trimmed read pairs; edit `config.sh` (paths, `activate_tools`).
2. `sbatch 01_strandedness_check.slurm` (run from this directory). Read `$WORKDIR/strand_check/result.txt`, then set `STRANDED=RF|FR|none` in `config.sh`.
3. `sbatch 02_hisat2_map.slurm` (array 1-8) -> `$WORKDIR/bams/`.
4. `sbatch --dependency=afterok:<jobid> 03_stringtie.slurm` (array 1-8) -> `$WORKDIR/stringtie/<sample>.gtf`.
5. `python3 04_extract_tss.py ref/Neonectria_faginata.gff3 tss_table.tsv stringtie/*.gtf`, then copy `tss_table.tsv` to `data/promoters/tss_table.tsv` here (small file; the BAMs stay on the cluster).
6. (optional, refine UTRs where StringTie kept the annotated start) `sbatch 06_split_strand.slurm` (array 1-8; strand-specific BAMs, needs STRANDED=RF/FR), then `sbatch --dependency=afterok:<jobid> 07_utr_coverage.slurm`, then `python3 08_merge_tss.py tss_table.tsv utr_coverage.tsv tss_table.merged.tsv` (prints StringTie-vs-coverage concordance). Copy `tss_table.merged.tsv` to `data/promoters/tss_table.tsv` (or pass its path to the classifier).

## Locally
`Rscript promoters/classify_snps_promoter.R` — uses `data/promoters/tss_table.tsv` if present, otherwise falls back to CDS coordinates (ATG as TSS, no UTR classes). Classes: coding (CDS span incl. introns), 5UTR, 3UTR, promoter_core (<=250 bp upstream of TSS), promoter_proximal (<=1 kb), 3prime_flank (<=500 bp downstream), distal_intergenic. Neighbour orientation is respected (convergent pairs have no promoter between them; divergent pairs flag `bidirectional`). Windows are parameters at the top of the script.
