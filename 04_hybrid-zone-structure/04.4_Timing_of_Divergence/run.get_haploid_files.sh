#!/bin/bash

name=MY_SAMPLE_NAME

for chrom in `seq 1 20`; do samtools mpileup -s -f ./panubis1_no_chr.fa -q30 -Q60 -r $chrom ./bams/${name}_MarkDuplicates.bam | /data/tunglab/tpv/Programs/pu2fa/pu2fa -c $chrom -C 40 > tmp.${name}.$chrom.haploid.fa; done

cat tmp.${name}*.haploid.fa > ./haploid_fa/$name.fa

rm tmp.${name}*.haploid.fa
