#!/bin/bash
path_genome=~/genomes/panubis1/Panubis_1.0.fa

# index will refer to the chromosome
index=${SLURM_ARRAY_TASK_ID}
## This will work as the chromosomes are indexed 1-20,X,Y,scafoldZZZ rather than "chr1", "chr2", etc.
chrom=$index

# load in relevant packages 
module load java/1.8.0_45-fasrc01
module load tabix
module load samtools
module load vcftools
export PATH=$PATH:~/Programs/cmake/bin/
module load gcc

# Get genotype calls
vcftools --gzvcf baboon_genotypes.raw.$chrom.vcf.gz  --keep 00_yel.list --max-missing 0 --max-alleles 2 --recode --out ./chrom_vcfs/02.yel.$chrom
vcftools --gzvcf baboon_genotypes.raw.$chrom.vcf.gz --keep 00_yellow_sources.list --max-missing 0 --max-alleles 2 --recode --out ./chrom_vcfs/02.yel_source.$chrom

vcftools --gzvcf baboon_genotypes.raw.$chrom.vcf.gz --keep 00_anu.list --max-missing 0 --max-alleles 2 --recode --out ./chrom_vcfs/02.anu.$chrom
vcftools --gzvcf baboon_genotypes.raw.$chrom.vcf.gz --keep 00_anubis_sources.list --max-missing 0 --max-alleles 2 --recode --out ./chrom_vcfs/02.anu_source.$chrom

# Merge genotypes using generate_gt
/data/tunglab/tpv/Programs/IBDmix/IBDmix-master/build/src/generate_gt -a ./chrom_vcfs/02.yel_source.$chrom.recode.vcf -m ./chrom_vcfs/02.yel.$chrom.recode.vcf -o ./chrom_vcfs/03.yellow.$chrom.gt
/data/tunglab/tpv/Programs/IBDmix/IBDmix-master/build/src/generate_gt -a ./chrom_vcfs/02.anu_source.$chrom.recode.vcf -m ./chrom_vcfs/02.anu.$chrom.recode.vcf -o ./chrom_vcfs/03.anubis.$chrom.gt

# run IBDmix, using each possible source in the "sources" file
for indiv in `cat 00_yellow_sources.list`; do ~/Programs/IBDmix/IBDmix-master/build/src/ibdmix -g ./chrom_vcfs/03.yellow.$chrom.gt --sample 00_yel.list --archaic $indiv --archaic-error 0.0025 -o IBDmix_by_chrom/yellow.relative_to_$indiv.$chrom.txt -t ; done
for indiv in `cat 00_anubis_sources.list`; do ~/Programs/IBDmix/IBDmix-master/build/src/ibdmix -g ./chrom_vcfs/03.anubis.$chrom.gt --sample 00_anu.list --archaic $indiv --archaic-error 0.0025 -o IBDmix_by_chrom/anubis.relative_to_$indiv.$chrom.txt -t ; done

# these are default parameters except the archaic is decreased to match the modern error rate. This is true because all of our genotypes are from high coverage modern samples, while the original use of this program was using Neanderthal genotype calls (hence the term archaic for the source individual). 
