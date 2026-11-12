#!/bin/bash
# One-time, on a login node with internet: conda env with StringTie + samtools (+ MEME for the later motif step).
# HISAT2 is taken from the cluster install, but it is also put in the env as a fallback.
set -euo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
conda create -y -n promoter -c conda-forge -c bioconda stringtie samtools hisat2 meme python=3.11
mkdir -p $HOME/Nf_promoters/{ref,logs,strand_check,bams,stringtie}
echo "Next: copy genome + gff3 into \$WORKDIR/ref, fill in \$WORKDIR/samples.tsv, edit config.sh"
