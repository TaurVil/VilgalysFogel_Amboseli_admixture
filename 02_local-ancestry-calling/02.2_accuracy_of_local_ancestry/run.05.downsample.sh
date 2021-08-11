#!/bin/bash

g=${SLURM_ARRAY_TASK_ID}

name=`head -$g 00names | tail -1`

module load samtools
samtools view -s 0.2 -b ./mapped_bams/$name.10x.bam > ./mapped_bams/$name.2x.bam

samtools view -s 0.1 -b ./mapped_bams/$name.10x.bam > ./mapped_bams/$name.1x.bam
