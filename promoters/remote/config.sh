#!/bin/bash
# Shared settings for the promoter / TSS pipeline. Edit the paths, then every script sources this file.
# Layout on the server (WORKDIR) is created by the scripts: ref/ hisat2_index/ strand_check/ bams/ stringtie/ logs/

export WORKDIR=$HOME/Nf_promoters                 # EDIT: working directory on the cluster
export GENOME=$WORKDIR/ref/Neonectria_faginata.scaffolds.fa   # copy from data/genotype_dat/
export GFF=$WORKDIR/ref/Neonectria_faginata.gff3              # copy from data/genotype_dat/
export SAMPLES=$WORKDIR/samples.tsv               # sample<TAB>R1<TAB>R2 (trimmed reads); see samples.template.tsv
export THREADS=8

# Library strandedness; set after running 01_strandedness_check.slurm.
#   RF   = dUTP / TruSeq stranded (read1 antisense to the gene)  -> hisat2 --rna-strandness RF, stringtie --rf
#   FR   = ligation-type stranded (read1 sense)                  -> hisat2 --rna-strandness FR, stringtie --fr
#   none = unstranded
export STRANDED=none

# Fungal introns are short; limits intron size for hisat2
export MAX_INTRON=3000

# EDIT to match your cluster: how to get hisat2, samtools, stringtie on PATH
activate_tools(){
    module purge
    module load anaconda/colsa                    # as in LFMM/sweepfinder2_scan.slurm; hisat2 lives here?
    # StringTie and samtools: conda env built by 00_env_setup.sh
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate promoter
}

hisat2_strand_flag(){ case "$STRANDED" in RF|FR) echo "--rna-strandness $STRANDED";; *) echo "";; esac; }
stringtie_strand_flag(){ case "$STRANDED" in RF) echo "--rf";; FR) echo "--fr";; *) echo "";; esac; }
