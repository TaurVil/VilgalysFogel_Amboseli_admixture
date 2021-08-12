#!/bin/bash

module load samtools

samp=SAMPLE_NAME

module load vcftools
vcftools --exclude-bed tracts_to_mask/to_remove.$samp.yes_intersect_50.bed --indv $samp --gzvcf ./refpanel.vcf.gz --recode --out indiv_vcfs/subvcf.$samp.$mask ; echo $samp

module load tabix
bgzip indiv_vcfs/subvcf.$samp.recode.vcf
tabix bgzip indiv_vcfs/subvcf.$samp.recode.vcf
