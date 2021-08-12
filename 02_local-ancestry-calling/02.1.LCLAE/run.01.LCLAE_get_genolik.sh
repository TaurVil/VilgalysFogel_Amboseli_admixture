#!/bin/bash

# run for each chromosome
index=${SLURM_ARRAY_TASK_ID}
chrom=$index

# only evaluate biallelic sites
module load vcftools
vcftools --max-alleles 2 --gzvcf amboseli_with_masked_refpanel.vcf.gz --chr $chrom --recode --out amboseli_with_masked_refpanel.max2.$chrom

module load tabix

module load bcftools
bcftools annotate -x ^INFO/DP,^FORMAT/GT,^FORMAT/AD,^FORMAT/DP,^FORMAT/GQ,^FORMAT/PL amboseli_with_masked_refpanel.max2.$chrom.recode.vcf | sed 's/|/\//g' > for_genolik.amboseli_with_masked_refpanel.max2.$chrom.vcf

# Replace XXX with the total number of individuals in the vcf (includes reference individuals -- e.g., unadmixed yellow and anubis -- and admixed individuals)
sed '/^#/d' for_genolik.amboseli_with_masked_refpanel.max2.$chrom.vcf | cut -f 1,2,10- | sed -e 's/:/ /g' | sed -e 's/\./999/g' | ./LCLAE/filtbaboon1b XXX > genolik.$chrom.genolik

# cut -f 1,2,10- --> get just chr, position, and genotypes from file
# sed -e 's/:/ /g' --> replace all ":" in genotypes with just spaces
# sed -e 's/\./999/g' --> replace all "." in genotypes with "999"
# sed -e 's/\// /g' --> remove all "/" between alleles in a genotype
# LCLAE/filtbaboon1b XXX <-- XXX specifies the TOTAL number of individuals in the file (includes reference individuals -- e.g., unadmixed yellow and anubis -- and admixed individuals to be tested)

mv genolik.$chrom.genolik genotype_likelihoods_maskedref/

rm for_genolik.amboseli_with_masked_refpanel.max2.$chrom.vcf
