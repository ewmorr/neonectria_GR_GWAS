#!/bin/bash
#SBATCH --job-name=utrcov
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=utrcov_%j.out
# Strand-specific pooled coverage upstream of each start codon / downstream of each stop codon -> UTR extension calls.
# Run after 06_split_strand.slurm has finished for all samples.
set -euo pipefail
source "$SLURM_SUBMIT_DIR/config.sh"
activate_tools
cd $WORKDIR
python3 $SLURM_SUBMIT_DIR/07_utr_from_coverage.py $GFF bams utr_coverage.tsv $(cut -f1 $SAMPLES | tr '\n' ' ')
