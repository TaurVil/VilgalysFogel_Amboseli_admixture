#!/bin/bash

# we ultimately need to convert our vcf into eigenstrat format in order to use ADMIXTOOLS/admixr
# first, convert vcf to plink format
module load vcftools
vcftools --gzvcf final.baboon.to.macam.n51.CHROM.vcf.gz --plink --out CHROM.n51

# second, update plink ped with correct population info for each sample
#--update-ids expects input with the following four fields:

#Old family ID
#Old within-family ID
#New family ID
#New within-family ID

# also need to convert our column 6 from weird 0 values to 1 (otherwise individuals are ignored)

module load plink

plink --file CHROM.n51 --update-ids updatepop.n51.txt --keep-allele-order --recode --out updated.n51.CHROM

plink --file updated.n51.CHROM --make-pheno names.pheno.txt '*' --allow-no-sex --keep-allele-order --recode --out updated2.n51.CHROM # names.pheno.txt is just the first two columns of updatepop.n51.txt

# second, convert plink format to eigenstrat format using the parameter file you set up
/data/tunglab/tpv/Programs/EIG-6.1.4/bin/convertf -p my.par.ped.eigenstrat.CHROM
