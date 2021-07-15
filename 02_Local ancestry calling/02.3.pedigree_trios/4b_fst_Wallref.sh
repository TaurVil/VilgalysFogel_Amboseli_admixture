#!/bin/bash

module load vcftools

# VCF containing all chromosomes of masked variants in yellow and anubis reference populations not yet merged with Amboseli
location=/data/tunglab/tpv/local_ancestry/unadmixed_individuals

vcftools --gzvcf $location/refpanel.vcf.gz --weir-fst-pop 00_anu_wall.list --weir-fst-pop 00_yel_wall.list.nodup --fst-window-size 35000 --fst-window-step 500 --out fst_unmasked_unmerged_Wallref_35kbwin_500bpstep
