# This directory contains an R script for evaluating whether genomic features predict the rate of change in anubis allele frequencies over time (Supplementary Methods 14.2).

## The R script, 1anubis_ancestry_change_over_time.R, requires two data sets, "amboseli.ancestry.250kwin.txt" and "amboseli.demographic.info.txt", which are available in the Duke Data Repository at XXXXXXXXX.
### amboseli.ancestry.250kwin.txt: includes anubis ancestry generated using masked SNPRC individuals for the reference panel (used in main text analyses) in 250 kb windows across the genome (excluding the 50 kb at the chromosome ends) for all Amboseli individuals.
### amboseli.demographic.info.txt: contains demographic data for Amboseli individuals with confirmed ids, including the starting and endings years they were present in the population. Note that other demographic data used in other scripts are also included in this data frame but which we won't use here (their sex, entry type (B = born into a study group, O = first observed in a new study group, or I = immigrating into a study group), genome_wide_anubis_ancestry)

## The R script also requires estimates of FST, mean recombination rate, and B-values. 
### FST: to get these estimates, use the script 4bFST_SNPRCref.sh in "VilgalysFogel_Amboseli_admixture/02_Local ancestry calling/02.3.pedigree_trios/" but substituting 250000 for window size --fst-window-size in place of 35000.
### Mean recombination rate and B-values: get these from ****XXXX TAURASSSSS TO_ANALYZE DATA FRAME***** which has already removed windows in which recombination rate estimates were greater than 100x larger than the median recombination rate.
