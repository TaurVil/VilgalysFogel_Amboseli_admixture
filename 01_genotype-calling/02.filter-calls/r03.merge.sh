#!/bin/bash

index=${SLURM_ARRAY_TASK_ID}

module load java/1.8.0_45-fasrc01
module load tabix
module load vcftools
module load bcftools


tabix ./chrom_vcfs/02.$index.amboseli.recode.vcf.gz
tabix chrom_vcfs/02.anu.$index.recode.vcf.gz
tabix chrom_vcfs/02.yel.$index.recode.vcf.gz

bcftools merge ./chrom_vcfs/02.$index.amboseli.recode.vcf.gz ./chrom_vcfs/02.anu.$index.recode.vcf.gz ./chrom_vcfs/02.yel.$index.recode.vcf.gz -O z -o chrom_vcfs/03.merged.$index.vcf.gz

vcftools --gzvcf ./chrom_vcfs/02.$index.amboseli.recode.vcf.gz --kept-sites --out 03.amboseli.$index
vcftools --gzvcf ./chrom_vcfs/02.anu.$index.recode.vcf.gz --kept-sites --out 03.reference.$index ## we already required the same sites show up in both yellow and anubis baboons

sort 03.amboseli.$index.kept.sites 03.reference.$index.kept.sites | uniq -d > 03.shared.$index.sites

vcftools --gzvcf chrom_vcfs/03.merged.$index.vcf.gz --positions 03.shared.$index.sites --recode --out 04.merged_shared.$index

mv 04.merged_shared.$index.recode.vcf 04.merged_shared.$index.vcf; bgzip 04.merged_shared.$index.vcf; tabix 04.merged_shared.$index.vcf.gz

rm chrom_vcfs/03.merged.$index.vcf.gz; rm 03.reference.$index.kept.sites ; rm 03.amboseli.$index.kept.sites

