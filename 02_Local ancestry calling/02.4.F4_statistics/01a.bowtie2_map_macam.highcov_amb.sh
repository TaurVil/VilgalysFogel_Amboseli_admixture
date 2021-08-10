#!/bin/bash

module load bowtie2
module load samtools

bowtie2 -p 16 -t -x MacaM_bowtie2/MacaM -1 *NAME*R1*.fastq.gz -2 *NAME*R2*.fastq.gz -S mapped.NAME.sam --no-unal

samtools view -bS mapped.NAME.sam > mapped.NAME.bam

rm mapped.NAME.sam
