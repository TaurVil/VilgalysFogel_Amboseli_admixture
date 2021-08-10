#!/bin/bash

# Need to convert vcfs into eigenstrat format in order to use ADMIXTOOLS/admixr
# To do so, convert vcf --> plink format --> eigenstrat format

module load vcftools
vcftools --gzvcf final.baboon.to.macam.n49.CHROM.rename.vcf.gz --mac 1 --plink --out CHROM.n49
# mac 1 = only include sites with a minor allele count of at least 1 (i.e., all samples are not all carrying the major allele) - doesn’t really matter because sites where individuals are completely for the major allele won’t be informative for f-statistics

# Second, update plink ped file with correct population info for each sample
#--update-ids expects input with the following four fields:

#Old family ID
#Old within-family ID
#New family ID
#New within-family ID

module load plink

plink --file CHROM.n49 --update-ids updatepop.n49.txt --keep-allele-order --recode --out updated.n49.CHROM

# Also, need to convert our column 6 from weird 0 values to 1 (otherwise individuals are ignored)
plink --file updated.n49.CHROM --make-pheno names.pheno.txt '*' --allow-no-sex --keep-allele-order --recode --out updated2.n49.CHROM # names.pheno.txt is just the first two columns of updatepop.n49.txt

# Finally, convert plink format to eigenstrat format using the convertf parameter
EIG-6.1.4/bin/convertf -p my.par.ped.eigenstrat.CHROM
