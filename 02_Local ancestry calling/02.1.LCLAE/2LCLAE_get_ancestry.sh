#!/bin/bash

# Run for each chromosome
index=${SLURM_ARRAY_TASK_ID}
chrom=$index

printf "`echo NUMBER` `echo $chrom` \n"

number=NUMBER

cat genotype_likelihoods_maskedref/genolik.$chrom.genolik | LCLAE/filtbaboon2c XXX ref_anubis.h  ref_yellow.h $number | LCLAE/geno_lik2 .2 35000 > $number.35kb.d2.$chrom.masked.SWref.txt

# LCLAE/filtbaboon2c XXX <-- XXX specifies the total number of individuals in file (e.g., if there are 100 individuals including reference and admixed individuals, put 100 in place of XXX)
# ref_anubis.h = reference population one (anubis baboons), ref_yellow.h = reference population two (yellow baboons) - LCLAE will estimate the number of alleles from reference population one (here, anubis) so 0 corresponds to 0 anubis alleles, 1 corresponds to 1 anubis allele, and 2 corresponds to 2 anubis alleles)
# LCLAE/geno_lik2 0.2 35000 <-- 0.2 specifies the minimum allele frequency difference (here, 20%) between the two reference populations, 35000 specifies the window size (here, 35 kb)
