#!/bin/bash

module load vcftools

# VCF containing all chromosomes of masked variants in yellow and anubis reference populations not yet merged with Amboseli individuals
location=/data/tunglab/tpv/panubis1_genotypes/masked_final

vcftools --gzvcf $location/masked.filtered.yes_intersect_50.vcf.gz --weir-fst-pop SNPRCanubis.list --weir-fst-pop SNPRCyellow.list --fst-window-size 35000 --fst-window-step 500 --out fst_masked_unmerged_SWref_35kbwin_500bpstep
