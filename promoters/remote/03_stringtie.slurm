#!/bin/bash
#SBATCH --job-name=stringtie
#SBATCH --array=1-8
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=stringtie_%A_%a.out
# Reference-guided StringTie per sample. Using -G keeps the funannotate gene IDs (Q2T73_...) for matched
# transcripts, so extended 5'/3' UTRs are attached to the genes already used in the LFMM analysis.
set -euo pipefail
source "$SLURM_SUBMIT_DIR/config.sh"
activate_tools
cd $WORKDIR; mkdir -p stringtie
read S R1 R2 < <(sed -n "${SLURM_ARRAY_TASK_ID}p" $SAMPLES)
stringtie bams/${S}.bam -p 4 -G $GFF $(stringtie_strand_flag) -l ${S} -o stringtie/${S}.gtf -A stringtie/${S}.abund.tsv
