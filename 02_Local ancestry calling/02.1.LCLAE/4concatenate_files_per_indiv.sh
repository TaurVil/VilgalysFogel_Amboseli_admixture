#!/bin/bash

# Concatenate all chromosome files per individual into a single file per individual
# h is the individual number
# for example, if you want to work on the first 10 individuals in the vcf, run:
for h in `seq 1 10`; do cat $h.35kb.d2.* > $h.35kb.d2.masked.SNPRCref.txt; done
