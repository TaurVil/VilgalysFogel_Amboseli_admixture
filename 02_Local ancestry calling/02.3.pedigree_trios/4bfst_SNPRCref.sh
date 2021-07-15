#!/bin/bash

module load vcftools

# VCF containing all chromosomes of masked variants in yellow and anubis reference populations not yet merged with Amboseli
location=/data/tunglab/tpv/panubis1_genotypes/masked_final

vcftools --gzvcf $location/masked.filtered.yes_intersect_50.vcf.gz --weir-fst-pop 00_anu_SW.list --weir-fst-pop 00_yel_SW.list --fst-window-size 35000 --fst-window-step 500 --out fst_masked_unmerged_SWref_35kbwin_500bpstep
