#!/bin/bash

module load bowtie2; module load samtools; module load tabix
g=${SLURM_ARRAY_TASK_ID}
name=`head -$g 00names | tail -1`

zcat ./sim_reads/$name.*fq.gz | bgzip > ./sim_reads/$name.fq.gz
tabix ./sim_reads/$name.fq.gz
#rm ./sim_reads/$name.chr*fq.gz

bowtie2 -p 16 -t -x ./reduced_genome -U ./sim_reads/$name.fq.gz --no-unal | samtools view -bS > mapped_bams/$name.10x.bam
# -U because this is SE data


