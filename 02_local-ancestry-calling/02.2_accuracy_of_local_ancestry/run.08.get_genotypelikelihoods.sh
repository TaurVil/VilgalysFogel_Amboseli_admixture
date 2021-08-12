#!/bin/bash

index=${SLURM_ARRAY_TASK_ID}
coverage=`head -$index 01_sets.txt | tail -1`
path_genome=./Reduced_Genome.fa

module load java; module load samtools; module load python
module load htslib
module load GATK
module load tabix

## Merge test samples with reference panel genotypes
java -jar ~/Programs/GenomeAnalysisTK.jar -T CombineVariants -R $path_genome --variant ./filt.$coverage.recode.vcf.gz --variant ../refpanel.vcf.gz -o merged.$coverage.vcf.gz -genotypeMergeOptions UNIQUIFY -L 01_targetted_chroms.bed
java -jar ~/Programs/GenomeAnalysisTK.jar -T SelectVariants -R $path_genome -V merged.$coverage.vcf.gz -select 'set == "Intersection"' -o CommonCalls.$coverage.vcf
vcftools --vcf CommonCalls.$coverage.vcf --max-alleles 2 --recode --out CommonCalls.biallelic.$coverage

## Extract genotype likelihoods
# LCLAE expects a certain set of information to properly split out the likelihoods 
bcftools annotate -x ^INFO/DP,^FORMAT/GT,^FORMAT/AD,^FORMAT/DP,^FORMAT/GQ,^FORMAT/PL CommonCalls.biallelic.$coverage.recode.vcf | sed 's/|/\//g' > merged.$coverage.forgenolik.vcf

sed '/^#/d' merged.$coverage.forgenolik.vcf | cut -f 1,2,10- | sed -e 's/:/ /g' | sed -e 's/\./999/g' | /data/tunglab/tpv/Programs/LCLAE/filtbaboon1b 56 > genolik.$coverage.genolik
## 70 is the total number of samples between reference individuals and individuals to call  (25 samples, 24 anubis, 7 yellow)
