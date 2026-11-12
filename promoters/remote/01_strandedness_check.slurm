#!/bin/bash
#SBATCH --job-name=strandchk
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=%x_%j.out
# Determine library strandedness: map 2M read pairs of the first sample UNSTRANDED, then ask whether read1
# lands on the same or opposite strand as the gene it overlaps (script check_strand.py).
set -euo pipefail
source "$SLURM_SUBMIT_DIR/config.sh"
activate_tools
cd $WORKDIR; mkdir -p hisat2_index strand_check
if [ ! -e hisat2_index/Nf.1.ht2 ]; then hisat2-build -p $THREADS $GENOME hisat2_index/Nf > hisat2_index/build.log 2>&1; fi
read S R1 R2 < <(head -1 $SAMPLES)
# first 2M pairs (8M lines)
zcat -f $R1 | head -8000000 > strand_check/sub_1.fq
zcat -f $R2 | head -8000000 > strand_check/sub_2.fq
hisat2 -p $THREADS --max-intronlen $MAX_INTRON -x hisat2_index/Nf -1 strand_check/sub_1.fq -2 strand_check/sub_2.fq 2> strand_check/hisat2.log \
  | samtools view -@ 4 -b -o strand_check/sub.bam -
samtools view -f 2 -F 256 strand_check/sub.bam | python3 "$SLURM_SUBMIT_DIR/check_strand.py" $GFF | tee strand_check/result.txt
rm -f strand_check/sub_1.fq strand_check/sub_2.fq
