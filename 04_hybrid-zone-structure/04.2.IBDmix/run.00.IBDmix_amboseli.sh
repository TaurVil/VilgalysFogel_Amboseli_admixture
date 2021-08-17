#!/bin/bash
path_genome=~/genomes/panubis1/Panubis_1.0.fa

# index will refer to the chromosome
index=${SLURM_ARRAY_TASK_ID}
## This will work as the chromosomes are indexed 1-20,X,Y,scafoldZZZ rather than "chr1", "chr2", etc.
chrom=$index

module load java/1.8.0_45-fasrc01
module load tabix
module load samtools
module load vcftools
export PATH=$PATH:~/Programs/cmake/bin/
module load gcc

# Get genotype calls
vcftools --gzvcf baboon_genotypes.raw.$chrom.vcf.gz --keep 00_amboseli.list --max-missing 0 --max-alleles 2 --recode --out ./chrom_vcfs/02.ambo.$chrom
vcftools --gzvcf baboon_genotypes.raw.$chrom.vcf.gz --keep 00_amboseli_sources.txt --max-missing 0 --max-alleles 2 --recode --out ./chrom_vcfs/02.source.$chrom

# Merge genotypes using generate_gt
~/Programs/IBDmix/IBDmix-master/build/src/generate_gt -a ./chrom_vcfs/02.source.$chrom.recode.vcf -m ./chrom_vcfs/02.ambo.$chrom.recode.vcf -o ./chrom_vcfs/03.$chrom.gt

# run IBDmix
for indiv in `cat 00_amboseli_sources.txt`; do ~/Programs/IBDmix/IBDmix-master/build/src/ibdmix -g ./chrom_vcfs/03.$chrom.gt --sample 00_amboseli.list --archaic $indiv --archaic-error 0.0025 -o IBDmix_by_chrom/amboseli.relative_to_$indiv.$chrom.txt -t ; done

