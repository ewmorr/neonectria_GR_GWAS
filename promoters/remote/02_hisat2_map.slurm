#!/bin/bash
#SBATCH --job-name=hisat2
#SBATCH --array=1-8
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=hisat2_%A_%a.out
# Map each sample (one array task per line of samples.tsv) and write a coordinate-sorted, indexed BAM.
# Set STRANDED in config.sh first (from 01_strandedness_check).
set -euo pipefail
source "$SLURM_SUBMIT_DIR/config.sh"
activate_tools
cd $WORKDIR; mkdir -p bams logs
[ -e hisat2_index/Nf.1.ht2 ] || hisat2-build -p $THREADS $GENOME hisat2_index/Nf > hisat2_index/build.log 2>&1
read S R1 R2 < <(sed -n "${SLURM_ARRAY_TASK_ID}p" $SAMPLES)
echo "sample $S  stranded=$STRANDED"
for f in "$R1" "$R2" hisat2_index/Nf.1.ht2; do [ -s "$f" ] || { echo "ERROR: missing or empty: $f"; exit 1; }; done
# On any failure show why: hisat2's stderr goes to a log file
trap 'echo "FAILED; tail of logs/${S}.hisat2.log:"; tail -n 20 logs/${S}.hisat2.log' ERR
hisat2 -p $THREADS --dta --max-intronlen $MAX_INTRON $(hisat2_strand_flag) -x hisat2_index/Nf -1 $R1 -2 $R2 \
    --summary-file logs/${S}.hisat2.summary 2> logs/${S}.hisat2.log \
  | samtools sort -@ 4 -m 2G -o bams/${S}.bam -
samtools index bams/${S}.bam
echo "done $S"; cat logs/${S}.hisat2.summary
