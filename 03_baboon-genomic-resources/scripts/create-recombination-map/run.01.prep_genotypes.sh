#!/bin/bash

# Filtering recommendations from Pfeifer 2020 Mol Biol Evol (vervet monkey genetic map)

## get SNPRC anubis from full list of baboon genome calls 
vcftools --gzvcf ./baboon1k_v1_snpEff_chr$chrom'.vcf.gz' --keep n24.anubis.SNPRCfounder.list --max-alleles 2 --remove-indels --recode --out ./tmp.1.$chrom.SNPRCanubis --recode-INFO-all
# vcftools filtering pass 1: keep biallelic variants called in all individuals 
# excluded sites with missing data to avoid errors and biases resulting from computational imputation
vcftools --gzvcf ./tmp.1.$chrom.SNPRCanubis.recode.vcf --max-missing 1 --minQ 30 --max-alleles 2 --remove-indels --maf 0.001 --recode --out ./tmp.2.$chrom.SNPRCanubis --recode-INFO-all

## GATK hard filtering 
## Version of the genome is one without the "chr" prefix for chromosome names, the same as was used in Section 01.2
java -jar ~/Programs/GenomeAnalysisTK.jar -T VariantFiltration -R ./Panubis1.nochromname.fa -V ./tmp.2.$chrom.SNPRCanubis.recode.vcf -filterName "FS" --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" -cluster 3 -window 10 -o ./tmp.3.$chrom.SNPRCanubis.vcf.gz

rm tmp.2.$chrom.SNPRCanubis*; rm tmp.1.$chrom.SNPRCanubis*


## Identify and remove singletons
## add in a minor allele frequency filter of 5%
vcftools --gzvcf ./tmp.3.$chrom.SNPRCanubis.vcf.gz --singletons --out ./SNPRCanubis.$chrom
vcftools --gzvcf ./tmp.3.$chrom.SNPRCanubis.vcf.gz --exclude-positions ./SNPRCanubis.$chrom.singletons --mag 0.05 --recode --out SNPRCanubis.$chrom --recode-INFO-all --remove-filtered-all
rm SNPRCanubis.*singletons 


## remove SNPs with excess heterozygosity
## Specifically, a P-value for Hardy–Weinberg Equilibrium was calculated using the “—hardy” option in VCFtools v.0.1.13 (Danecek et al. 2011), and SNPs with P < 0.01 removed.
vcftools --hardy --gzvcf SNPRCanubis.$chrom.recode.vcf --out SNPRCanubis.$chrom
module load R; R 
library(data.table); fread("SNPRCanubis.$chrom.hwe") -> hwe; subset(hwe, hwe$P_HWE < 0.01) -> hwe; dim(hwe); write.table(hwe[,1:2], "anubisSNPRC_hwe_to_remove.$chrom.sites", row.names=F, col.names=F, sep="\t", quote=F); quit(save="no")
vcftools --gzvcf SNPRCanubis.$chrom.recode.vcf --recode-INFO-all --remove-filtered-all --exclude-positions anubisSNPRC_hwe_to_remove.$chrom.sites --recode --out SNPRCanubis.$chrom.hwe 

bgzip SNPRCanubis.$chrom.hwe.recode.vcf ; tabix SNPRCanubis.$chrom.hwe.recode.vcf.gz


## Pfeifer 2020 only used SNPs that could be reciprocally lifted over with the human genome. We'll not include this filter, assuming the quality of the baboon genome is sufficient. 
## Remove fixed alleles: done earlier with the maf frequency
