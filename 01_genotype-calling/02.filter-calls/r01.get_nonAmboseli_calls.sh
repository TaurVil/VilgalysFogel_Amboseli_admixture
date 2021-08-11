#!/bin/bash

## Extract unadmixed yellow and anubis genotypes

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


# Get genotypes called in either reference panel, then remove invariant sites and prepare to remove singletons
vcftools --gzvcf ./baboon1k_v1_snpEff_chr$chrom'.vcf.gz' --keep 00_anu.list --keep 00_yel.list --max-alleles 2 --remove-indels --recode --out ./chrom_vcfs/01.raw_biallelic_snps.$chrom.unadmixed --recode-INFO-all
vcftools --gzvcf ./chrom_vcfs/01.raw_biallelic_snps.$chrom.unadmixed.recode.vcf --max-missing 0.5 --minQ 30 --max-alleles 2 --remove-indels --maf 0.001 --recode --out ./chrom_vcfs/01.$chrom.missing_minQ_biallelic.unadmixed --recode-INFO-all
vcftools --gzvcf ./chrom_vcfs/01.$chrom.missing_minQ_biallelic.unadmixed.recode.vcf --singletons --out ./chrom_vcfs/01.$chrom.unadmixed

## GATK hard filtering criteria. Intentionally liberal.
## https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants
## Panubis1.nochromname.fa is a version of the baboon genome where chromosomes are index 1-20 rather than "chr1", "chr2", etc. 
java -jar ~/Programs/GenomeAnalysisTK.jar -T VariantFiltration -R ./Panubis1.nochromname.fa -V ./chrom_vcfs/01.$chrom.missing_minQ_biallelic.unadmixed.recode.vcf -filterName "FS" --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" -o ./chrom_vcfs/01b.gatk_filtered.unadmixed.$chrom.vcf.gz

# Get sitest that should be removed due to insufficient species-specific coverage
vcftools --gzvcf ./chrom_vcfs/01.$chrom.missing_minQ_biallelic.unadmixed.recode.vcf --keep 00_anu.list --max-missing 0.5 --remove-indels --removed-sites --out ./chrom_vcfs/01.anu.$chrom
vcftools --gzvcf ./chrom_vcfs/01.$chrom.missing_minQ_biallelic.unadmixed.recode.vcf --keep 00_yel.list --max-missing 0.5 --remove-indels --removed-sites --out ./chrom_vcfs/01.yel.$chrom
# Merge all sites to remove (singletons and low species coverage)
cat ./chrom_vcfs/01.$chrom.unadmixed.singletons ./chrom_vcfs/01.anu.$chrom.removed.sites ./chrom_vcfs/01.yel.$chrom.removed.sites > ./chrom_vcfs/01.$chrom.to_remove.txt


## output separate vcf, retaining only non-singleton sites with adequate coverage
vcftools --gzvcf ./chrom_vcfs/01b.gatk_filtered.unadmixed.$chrom.vcf.gz --exclude-positions ./chrom_vcfs/01.$chrom.to_remove.txt --keep 00_anu.list --recode --out ./chrom_vcfs/02.anu.$chrom --recode-INFO-all --remove-filtered-all
vcftools --gzvcf ./chrom_vcfs/01b.gatk_filtered.unadmixed.$chrom.vcf.gz --exclude-positions ./chrom_vcfs/01.$chrom.to_remove.txt --keep 00_yel.list --recode --out ./chrom_vcfs/02.yel.$chrom --recode-INFO-all --remove-filtered-all

# clean up and compress all output files 
bgzip ./chrom_vcfs/01.$chrom.missing_minQ_biallelic.unadmixed.recode.vcf
bgzip ./chrom_vcfs/02.yel.$chrom.recode.vcf
bgzip ./chrom_vcfs/02.anu.$chrom.recode.vcf
bgzip ./chrom_vcfs/01.raw_biallelic_snps.$chrom.unadmixed.recode.vcf

