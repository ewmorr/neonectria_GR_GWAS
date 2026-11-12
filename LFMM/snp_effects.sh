#!/bin/bash
# Predict coding consequences (synonymous / missense / ...) for every SNP in the LFMM VCF using bcftools csq.
# Run from repo root, with bcftools on PATH (conda activate bcftools).
# Output: data/LFMM/snp_effects/all_snps.csq.txt  (CHROM POS REF ALT BCSQ; contig names WITHOUT _pilon)
set -euo pipefail
G=data/genotype_dat
VCF=$G/FINAL_snp.mac_ge2.biallele.retainNG152.gwas_analyses.vcf.gz
OUT=data/LFMM/snp_effects
mkdir -p $OUT

# 1. csq needs Ensembl-style GFF3: gene:/transcript: prefixed IDs + biotype. Protein-coding (mRNA) models only.
awk -F'\t' '$3=="mRNA"{match($9,/ID=[^;]+/);print substr($9,RSTART+3,RLENGTH-3)}' $G/Neonectria_faginata.gff3 > $OUT/tx.txt
awk -F'\t' -v OFS='\t' -v TX=$OUT/tx.txt 'BEGIN{while((getline l<TX)>0)tx[l]=1}
/^#/{next}
{ id="";par="";
  if(match($9,/ID=[^;]+/))id=substr($9,RSTART+3,RLENGTH-3);
  if(match($9,/Parent=[^;]+/))par=substr($9,RSTART+7,RLENGTH-7);
  if($3=="mRNA"){ $9="ID=transcript:"id";Parent=gene:"par";biotype=protein_coding"; print }
  else if($3=="gene"){ if((id"-T1") in tx){ $9="ID=gene:"id";Name="id";biotype=protein_coding"; print } }
  else if($3=="exon"||$3=="CDS"){ if(par in tx){ $9="Parent=transcript:"par; print } }
}' $G/Neonectria_faginata.gff3 > $OUT/csq.gff3

# 2. sites-only VCF with contig names matching the fasta/gff (strip _pilon)
bcftools view -G $VCF -Ou | bcftools query -f '%CHROM\n' | sort -u | awk '{print $1"\t"substr($1,1,length($1)-6)}' > $OUT/rename.txt
bcftools view -G $VCF -Ou | bcftools annotate --rename-chrs $OUT/rename.txt -Oz -o $OUT/all_snps.sites.vcf.gz

# 3. -l/--local-csq: treat each SNP independently (no haplotype/phase logic; haploid, sites-only)
bcftools csq -l -f $G/Neonectria_faginata.scaffolds.fa -g $OUT/csq.gff3 --force -Ov $OUT/all_snps.sites.vcf.gz \
  | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/BCSQ\n' > $OUT/all_snps.csq.txt
wc -l $OUT/all_snps.csq.txt
