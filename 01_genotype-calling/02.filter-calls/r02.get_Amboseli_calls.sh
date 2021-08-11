#!/bin/bash

## Extract amboseli genotypes

## chromsomes are indexed 1-20, X,Y rather than "chr1", "chr2", etc. 
index=${SLURM_ARRAY_TASK_ID}
chrom=$index

## load relevant packages
module load java/1.8.0_45-fasrc01
module load tabix
module load samtools
module load vcftools
export PATH=$PATH:~/Programs/cmake/bin/
module load gcc

## extract Amboseli samples from full vcf file
## get just biallelic SNPs, not indels 
vcftools --gzvcf ./baboon1k_v1_snpEff_chr$chrom'.vcf.gz' --keep 00_amboseli.list --max-alleles 2 --remove-indels --recode --out ./chrom_vcfs/01.raw_biallelic_snps.$chrom.amboseli --recode-INFO-all

## GATK hard filtering criteria. Intentionally liberal.
## https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants
## Panubis1.nochromname.fa is a version of the baboon genome where chromosomes are index 1-20 rather than "chr1", "chr2", etc. 
java -jar ~/Programs/GenomeAnalysisTK.jar -T VariantFiltration -R ./Panubis1.nochromname.fa -V ./chrom_vcfs/01.raw_biallelic_snps.$chrom.amboseli.recode.vcf -filterName "FS" --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" -o ./chrom_vcfs/01b.gatk_filtered.amboseli.$chrom.vcf.gz

## Filter for sites called in 20% of Amboseli individuals, and variable within Amboseli (the maf threshold)
vcftools --gzvcf ./chrom_vcfs/01b.gatk_filtered.amboseli.$chrom.vcf.gz --max-missing 0.2 --minQ 30 --max-alleles 2 --remove-indels --maf 0.001 --recode --out ./chrom_vcfs/02.$chrom.amboseli --recode-INFO-all --remove-filtered-all

bgzip ./chrom_vcfs/01.raw_biallelic_snps.$chrom.amboseli.recode.vcf
bgzip ./chrom_vcfs/02.$chrom.amboseli.recode.vcf

