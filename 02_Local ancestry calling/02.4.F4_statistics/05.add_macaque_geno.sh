#!/bin/bash

module load bcftools

# Remove vcf header
bcftools view --no-header baboon.to.macam.n48.CHROMOSOME.filtused.recode.vcf >> noheader.baboon.to.macam.n48.CHROMOSOME.filtused.recode.vcf
# Add a new sample for all sites where the genotype is homozygous reference i.e., 0/0 (other genotype fields - e.g., AD:DP:GQ:PL - were grabbed from a sample with the homozygous reference genotype but these fields will not be important later)
sed 's/$/\t0\/0:49,0:49:99:0,120,1800/g' noheader.baboon.to.macam.n48.CHROMOSOME.filtused.recode.vcf >> tmp.CHROMOSOME.vcf

# Add macaque as the name of our newly added sample
sed '/^#CHROM/ {s/$/\tmacaque/}' baboon.to.macam.n48.CHROMOSOME.filtused.recode.vcf >> tmp2.CHROMOSOME.vcf

# Bgzip and index our vcf with the updated sample name in the header
module load tabix

bgzip tmp2.CHROMOSOME.vcf
tabix tmp2.CHROMOSOME.vcf.gz

# Only output the updated vcf header
bcftools view -h tmp2.CHROMOSOME.vcf.gz >> header.CHROMOSOME

# Add new header to our new vcf
cat header.CHROMOSOME tmp.CHROMOSOME.vcf >> final.baboon.to.macam.n49.CHROMOSOME.vcf

# Bgzip and index our new vcf (n49 since we added one more sample - the macaque)
bgzip final.baboon.to.macam.n49.CHROMOSOME.vcf  
tabix final.baboon.to.macam.n49.CHROMOSOME.vcf.gz

rm tmp*CHROMOSOME*vcf*
rm *head*CHROMOSOME*
