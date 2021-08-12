#!/bin/bash
module load tabix; module load samtools; module load vcftools; module load bcftools

mask=VERSION_NAME

## merge vcf files for each individual
bcftools merge ./indiv_vcfs/subvcf.*.$mask.recode.vcf.gz -O z -o ./refpanel_vcfs/masked.$mask.vcf.gz; echo $mask

## identify singletons to remove
vcftools --gzvcf ./refpanel_vcfs/masked.$mask.vcf.gz --singletons --out ./01.$mask.unadmixed
## identify sites that no longer fit coverage criteria (i.e. too much missing data)
vcftools --gzvcf ./refpanel_vcfs/masked.$mask.vcf.gz --keep ../../panubis1_genotypes/00_anu.list --max-missing-count 62 --remove-indels --removed-sites --out ./01.anu.$mask
vcftools --gzvcf ./refpanel_vcfs/masked.$mask.vcf.gz --keep ../../panubis1_genotypes/00_yel.list --max-missing-count 24 --remove-indels --removed-sites --out ./01.yel.$mask
## combine list of sites to remove
cat ./01.$mask.unadmixed.singletons ./01.anu.$mask.removed.sites ./01.yel.$mask.removed.sites > ./01.$mask.to_remove.txt
rm ./01.$mask.unadmixed.singletons ./01.anu.$mask.removed.sites ./01.yel.$mask.removed.sites

## get filtered vcf file 
vcftools --gzvcf ./refpanel_vcfs/masked.$mask.vcf.gz --exclude-positions ./01.$mask.to_remove.txt --recode --out ./refpanel_vcfs/masked.filtered.$mask --recode-INFO-all --remove-filtered-all

## rename and cleanup
mv ./refpanel_vcfs/masked.filtered.$mask.recode.vcf ./masked_yellow_and_anubis.vcf; bgzip ./masked_yellow_and_anubis.vcf; tabix ./masked_yellow_and_anubis.vcf.gz
rm ./01.$mask.to_remove.txt; rm -r ./refpanel_vcfs/; rm -r ./indiv_vcfs/