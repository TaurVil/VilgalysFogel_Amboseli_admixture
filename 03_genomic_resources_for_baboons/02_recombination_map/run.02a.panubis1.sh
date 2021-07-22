#!/bin/bash

module load vcftools; module load tabix; module load samtools

module load gsl; module load gcc/4.9.3-fasrc01; module load boost/1.54.0-fasrc01; module load python/2.7.11-fasrc01

cat individual_haplotype_fastas/*.CHROM.*.fa > chrCHROM.fa
/data/tunglab/tpv/Programs/LDhelmet_v1.10/ldhelmet find_confs --num_threads 24 -w 50 -o output/outputCHROM.conf chrCHROM.fa
time /data/tunglab/tpv/Programs/LDhelmet_v1.10/ldhelmet table_gen --num_threads 48 -t 0.0016 -r 0.0 0.1 2.0 0.5 10.0 5.0 100.0 -c output/outputCHROM.conf -o output/outputCHROM.lk
time /data/tunglab/tpv/Programs/LDhelmet_v1.10/ldhelmet pade --num_threads 48 -t 0.0016 -x 11 --defect_threshold 40 -c output/outputCHROM.conf -o output/outputCHROM.pade

time /data/tunglab/tpv/Programs/LDhelmet_v1.10/ldhelmet rjmcmc --num_threads 48 -l output/outputCHROM.lk -p output/outputCHROM.pade -s chrCHROM.fa -b 5 --burn_in 10000 -n 100000 -o output/outputCHROM.post
time /data/tunglab/tpv/Programs/LDhelmet_v1.10/ldhelmet post_to_text -m -p 0.025 -p 0.50 -p 0.975 -o output/output.CHROM.txt output/outputCHROM.post

