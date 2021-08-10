#!/bin/bash

module load bowtie2
module load samtools

# for fastq files where paired-end data are labeled as _1 and _2
bowtie2 -p 16 -t -x MacaM_bowtie2/MacaM -1 NAME*_1*.fastq.gz -2 NAME*_2*.fastq.gz -S mapped.NAME.sam --no-unal

samtools view -bS mapped.NAME.sam > mapped.NAME.bam

rm mapped.NAME.sam
