#!/bin/bash

chrom=CHROMOSOME

printf "`echo NUMBER` `echo $chrom` \n"

number=NUMBER

cat ./genotype_likelihoods_maskedref/genolik.$chrom.genolik | LCLAE/filtbaboon2c 508 ref_anubis.h ref_yellow.h $number | LCLAE/geno_lik2 .2 35000 > $number.35kb.d2.$chrom.maskedSWref.txt
