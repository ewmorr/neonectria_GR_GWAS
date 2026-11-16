#!/bin/bash
#SBATCH --job-name=splitstrand
#SBATCH --array=1-8
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=splitstrand_%A_%a.out
# Split each sample's BAM by TRANSCRIPT strand (primary alignments only) for strand-specific coverage.
#   RF (dUTP):  transcript on + strand  <=> read1 maps reverse  or  read2 maps forward
#   FR:         transcript on + strand  <=> read1 maps forward  or  read2 maps reverse
set -euo pipefail
source "$SLURM_SUBMIT_DIR/config.sh"
activate_tools
cd $WORKDIR
read S R1 R2 < <(sed -n "${SLURM_ARRAY_TASK_ID}p" $SAMPLES)
case "$STRANDED" in
  RF) PLUS_A="-f 80  -F 260"; PLUS_B="-f 128 -F 276"; MINUS_A="-f 64 -F 276"; MINUS_B="-f 144 -F 260";;
  FR) PLUS_A="-f 64  -F 276"; PLUS_B="-f 144 -F 260"; MINUS_A="-f 80 -F 260"; MINUS_B="-f 128 -F 276";;
  *)  echo "ERROR: strand-specific coverage needs STRANDED=RF or FR (got $STRANDED)"; exit 1;;
esac
T=$(mktemp -d); trap 'rm -rf $T' EXIT
samtools view -@ 4 -b $PLUS_A  bams/${S}.bam -o $T/pa.bam;  samtools view -@ 4 -b $PLUS_B  bams/${S}.bam -o $T/pb.bam
samtools merge -f -@ 4 bams/${S}.plus.bam  $T/pa.bam $T/pb.bam
samtools view -@ 4 -b $MINUS_A bams/${S}.bam -o $T/ma.bam;  samtools view -@ 4 -b $MINUS_B bams/${S}.bam -o $T/mb.bam
samtools merge -f -@ 4 bams/${S}.minus.bam $T/ma.bam $T/mb.bam
samtools index bams/${S}.plus.bam; samtools index bams/${S}.minus.bam
echo "$S plus reads: $(samtools view -c bams/${S}.plus.bam)  minus reads: $(samtools view -c bams/${S}.minus.bam)"
