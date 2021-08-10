#!/bin/bash

module load samtools

samtools merge mapped.SAMPLE1.bam bams/sort.SAMPLE2.bam bams/sort.SAMPLE3.bam bams/sort.SAMPLE4.bam bams/sort.SAMPLE5.bam
