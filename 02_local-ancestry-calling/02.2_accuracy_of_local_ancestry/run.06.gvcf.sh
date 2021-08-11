#!/bin/bash

module load java; module load samtools

g=${SLURM_ARRAY_TASK_ID}

f=`head -$g 00_bams.txt | tail -1`;

path_genome=./Reduced_Genome.fa

samtools addreplacerg -o tmp.$f.bam -r ID:$f -r PL:Illumina -r SM:$f mapped_bams/$f.bam
mv tmp.$f.bam mapped_bams/${f}.bam

samtools sort -o tmp.$f.bam mapped_bams/${f}.bam; mv tmp.$f.bam mapped_bams/${f}.bam

rm mapped_bams/${f}.bam.bai

samtools index mapped_bams/${f}.bam

module load GATK

GenomeAnalysisTK.sh HaplotypeCaller -ERC GVCF -I mapped_bams/$f.bam  -R $path_genome  -O gVCF/$f.g.vcf.gz -mbq 20

