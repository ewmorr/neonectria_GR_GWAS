conda activate bcftools

cd ~/repo/neonectria_GR_GWAS/data/

# filter the biallele table to just the GWAS samples
bcftools view -S sample_metadata/sequence_IDs.Nf.txt genotype_dat/FINAL_snp.biallele.retainNG152.vcf.gz -Oz -o genotype_dat/FINAL_snp.biallele.GWAS_samples.vcf.gz

######
cd genotype_dat
# convert to alelle count table
vcftools --gzvcf FINAL_snp.biallele.GWAS_samples.vcf.gz --counts2 --out FINAL_snp.biallele.GWAS_samples
# kept 969175 out of a possible 969175 Sites
# 
# for seepfinder we need to filter out sites with count == 0 (i.e., monomorphic). We also need to split by chromosome to retain the correct position argument
head FINAL_snp.biallele.GWAS_samples.frq.count
tail -n+2 FINAL_snp.biallele.GWAS_samples.frq.count | awk -v OFS="\t" '{print $2,$6,$4,"1"}' > FINAL_snp.biallele.GWAS_samples.sf2
echo -e 'position\tx\tn\tfolded' | cat - FINAL_snp.biallele.GWAS_samples.sf2 > temp && mv temp FINAL_snp.biallele.GWAS_samples.sf2
head FINAL_snp.biallele.GWAS_samples.sf2

head FINAL_snp.biallele.GWAS_samples.frq.count | grep -v '\s0$' 
gunzip -c FINAL_snp.biallele.GWAS_samples.vcf.gz | less -S

grep -v '\s0$' FINAL_snp.biallele.GWAS_samples.frq.count > FINAL_snp.biallele.GWAS_samples.nozeros.frq.count
tail -n+2 FINAL_snp.biallele.GWAS_samples.nozeros.frq.count | awk -v OFS="\t" '{print $2,$6,$4,"1"}' > FINAL_snp.biallele.GWAS_samples.sf2
echo -e 'position\tx\tn\tfolded' | cat - FINAL_snp.biallele.GWAS_samples.sf2 > temp && mv temp FINAL_snp.biallele.GWAS_samples.sf2

# for smoe reason the below doesn't work in RStudio IDE shell so run in nomral shell
while IFS= read -r line
do(
    echo $line
    grep $line FINAL_snp.biallele.GWAS_samples.frq.count | \
        grep -v '\s0$' | \
        tail -n+2 | \
        awk -v OFS="\t" '{print $2,$6,$4,"1"}' > $line.sf2
        echo -e 'position\tx\tn\tfolded' | cat - $line.sf2 > temp && mv temp $line.sf2
)
done < tigs_names.txt
mkdir sf2_files
mv *sf2 sf2_files


######
# on premise
# 

SweepFinder2 -f FINAL_snp.biallele.GWAS_samples.sf2 SpecFile

module load linuxbrew/colsa
while IFS= read -r line
do(
    echo $line
    SweepFinder2 –f sf2_files/FINAL_snp.biallele.GWAS_samples.$line.sf2 spect_files/$line.spct
    
)
done < tigs_names.txt
# 