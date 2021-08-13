#!/bin/bash

# process mapped samples using samtools and picard tools
# sort, add read groups, remove duplicates, and remove alignments with mapping quality less than 10
module load java
module load samtools

# Sort data
samtools sort -o sort.SAMPLE.bam ./bams/mapped.SAMPLE.bam

# Add read groups
java -jar picard-tools-1.137/picard.jar AddOrReplaceReadGroups VALIDATION_STRINGENCY=LENIENT INPUT=./sort.SAMPLE.bam OUTPUT=./rg.SAMPLE.bam RGPL=illumina RGLB=reseq RGPU=reseq RGSM=SAMPLE

samtools index rg.SAMPLE.bam

rm sort.SAMPLE.bam

# Remove duplicates
java -jar picard-tools-1.137/picard.jar MarkDuplicates I=rg.SAMPLE.bam O=nodup.SAMPLE.bam CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true M=nodup.SAMPLE.metrics

rm rg.SAMPLE.bam

# Remove alignments with mapping quality less than 10
samtools view -q 10 -o nodup.mapq10.SAMPLE.bam nodup.SAMPLE.bam
samtools index nodup.mapq10.SAMPLE.bam

rm nodup.SAMPLE.bam
