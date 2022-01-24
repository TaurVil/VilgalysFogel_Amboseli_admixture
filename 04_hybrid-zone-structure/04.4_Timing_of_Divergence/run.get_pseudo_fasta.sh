#!/bin/bash

yelname=MY_YELLOW_NAME
anuname=MY_ANUBIS_NAME

module load python/2.7.6-fasrc01

python ../../Programs/hPSMC/psmcfa_from_2_fastas.py -b10 -m5 haploid_fa/$yelname.fa haploid_fa/$anuname.fa > ./$yelname.$anuname.psmcfa

mv ./$yelname.$anuname.psmcfa ./pseudo_fasta/

