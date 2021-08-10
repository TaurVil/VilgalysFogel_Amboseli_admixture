#!/bin/bash

# Need to convert vcfs into eigenstrat format in order to use ADMIXTOOLS/admixr
# To do so, convert vcf --> plink format --> eigenstrat format

# First, convert vcf to plink format
module load vcftools
vcftools --gzvcf final.baboon.to.macam.n48.CHROM.vcf.gz --plink --out CHROM.n51

# Second, update plink ped file with correct population info for each sample
#--update-ids expects input with the following four fields:

#Old family ID
#Old within-family ID
#New family ID
#New within-family ID

module load plink

plink --file CHROM.n51 --update-ids updatepop.n51.txt --keep-allele-order --recode --out updated.n51.CHROM

# Also, need to convert our column 6 from weird 0 values to 1 (otherwise individuals are ignored)
plink --file updated.n51.CHROM --make-pheno names.pheno.txt '*' --allow-no-sex --keep-allele-order --recode --out updated2.n51.CHROM # names.pheno.txt is just the first two columns of updatepop.n48.txt

# Finally, convert plink format to eigenstrat format using the convertf parameter
EIG-6.1.4/bin/convertf -p my.par.ped.eigenstrat.CHROM
