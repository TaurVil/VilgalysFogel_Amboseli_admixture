#!/bin/bash

module load vcftools

# VCF containing all chromosomes of masked variants in yellow and anubis reference populations not yet merged with Amboseli individuals
# this VCF file is available from Zenodo

vcftools --gzvcf masked_yellow_and_anubis.vcf.gz --weir-fst-pop Wallanubis.list --weir-fst-pop Wallyellow.list --fst-window-size 35000 --fst-window-step 500 --out fst_unmasked_unmerged_Wallref_35kbwin_500bpstep
