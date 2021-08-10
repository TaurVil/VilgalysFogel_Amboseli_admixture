#!/bin/bash

module load samtools

samtools sort -o sort.NAME.bam bams/mapped.NAME.bam
