#!/bin/bash

for f in `cat 00_vcf_sample_order.masked.list`; do awk 'BEGIN {OFS="\t"} {print $0,FILENAME}' $f.35kb.d2.masked.SNPRCref.txt > tmp.$f.txt; mv tmp.$f.txt $f.35kb.d2.masked.SNPRCref.txt; sed -i 's/.35kb.d2.masked.SWref.txt//g' $f.35kb.d2.masked.SWref.txt; done
