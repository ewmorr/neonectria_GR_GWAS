#!/bin/bash
#SBATCH --job-name=sf2
#SBATCH --output=sweep_array.out
##SBATCH --error=logs/%x_%A_%a.err
#SBATCH --array=0-23%8
#SBATCH --cpus-per-task=1
#SBATCH --exclude=node117,node118

module purge
module load linuxbrew/colsa

workdir=$HOME/Nf_sweepfinder

freq_files=(${workdir}/tig*.sf2)
freq=${freq_files[$SLURM_ARRAY_TASK_ID]}
tig=$(basename "$freq" _pilon.sf2)

echo "Starting $tig"

SweepFinder2 -lg 12500 $freq SpecFile out_${tig}.txt
