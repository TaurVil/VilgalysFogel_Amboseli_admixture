#!/bin/bash

module load java
module load samtools

g=${SLURM_ARRAY_TASK_ID}

f=`head -$g samples_final | tail -1`;

path_genome=MacaM_bowtie2/MacaM.fa

module load GATK

GenomeAnalysisTK.sh HaplotypeCaller -ERC GVCF -I bams/nodup.mapq10.$f.bam -R $path_genome -L CHROM -O gVCF/$f.CHROM.g.vcf.gz -mbq 20
#-mbq 20 = Minimum base quality required to consider a base for calling of 20
