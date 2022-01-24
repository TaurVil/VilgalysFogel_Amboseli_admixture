#!/bin/bash

yelname=MY_YELLOW_NAME
anuname=MY_ANUBIS_NAME

module load python/2.7.6-fasrc01

../../Programs/PSMC/psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o $yelname.$anuname.psmc pseudo_fasta/$yelname.$anuname.psmcfa

