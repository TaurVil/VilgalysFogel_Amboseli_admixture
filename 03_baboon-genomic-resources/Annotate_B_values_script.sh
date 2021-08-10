#!/bin/bash
#SBATCH --get-user-env

module load gsl
module load gcc

## Parameters the results of optimization in McVicker et al. 2009, fit to the density of variable sites between the human and chimpanzee genomes. 

## let's do exonic first
u=7.4e-8 # mutation rate 
a=2.5
b=-3
feature=exons

#sed -e s/SEDu/$u/g 00_panubis1_to_sed.conf | sed -e s/SEDa/$a/g | sed -e s/SEDb/$b/g | sed -e s/SEDFEATURE/$feature/g > tmp.$index.conf
#mkdir ./output.$u.$a.$b.$feature/
#../bkgd-master/calc_bkgd tmp.$index.conf
#rm tmp.$index.conf

## now time for non-exonic
u=8.4e-10
a=1
b=-5
feature=genes_plus_10kbprom

index=1401

sed -e s/SEDu/$u/g 00_panubis1_to_sed.conf | sed -e s/SEDa/$a/g | sed -e s/SEDb/$b/g | sed -e s/SEDFEATURE/$feature/g > tmp.$index.conf
mkdir ./output.$u.$a.$b.$feature/
../bkgd-master/calc_bkgd tmp.$index.conf
rm tmp.$index.conf
