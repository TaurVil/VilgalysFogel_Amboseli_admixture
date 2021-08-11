#!/bin/bash

path_genome=MacaM_bowtie2/MacaM.fa

module load java; module load samtools; module load python
module load htslib
module load GATK
module load tabix

ls gVCF/*.CHROM.g.vcf.gz > indivs.CHROM.list

GenomeAnalysisTK.sh CombineGVCFs  -R $path_genome -O baboon.to.macam.n94.CHROM.g.vcf.gz -L CHROM -V indivs.CHROM.list
GenomeAnalysisTK.sh GenotypeGVCFs -R $path_genome -V baboon.to.macam.n94.CHROM.g.vcf.gz -L CHROM -O baboon.to.macam.n94.CHROM.vcf.gz 

tabix -f baboon.to.macam.n94.CHROM.vcf.gz
module load java/1.8.0_45-fasrc01
module load vcftools

# GATK hard filtering for high quality variants, remove clusters of 3 or more variants that fell within a 10 bp window
java -jar /data/tunglab/tpv/Programs/GenomeAnalysisTK.jar -T VariantFiltration -R $path_genome -V baboon.to.macam.n94.CHROM.vcf.gz -filterName "FS" --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" -cluster 3 -window 10 -o baboon.to.macam.n94.CHROM.filt.vcf.gz

# Only keep the 48 individuals we will analyze for f4 stats
#vcftools --gzvcf baboon.to.macam.n94.CHROM.filt.vcf.gz --keep f4indivlist_n48  --recode --out tmp.n48.CHROM --recode-INFO-all

# Identify singletons
vcftools --vcf tmp.n48.CHROM.recode.vcf --singletons --out tmp.n48.CHROM

# Remove singletons, doubletons, indels, and sites not typed in all individuals; retain only biallelic sites 
vcftools --vcf tmp.n48.CHROM.recode.vcf --exclude-positions tmp.n48.CHROM.singletons --max-alleles 2 --max-missing 1 --remove-indels --recode --out baboon.to.macam.n48.CHROM.filtused --recode-INFO-all --remove-filtered-all

rm tmp.n48.CHROM.recode.vcf
