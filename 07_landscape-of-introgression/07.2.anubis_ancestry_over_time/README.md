## This directory contains an R script for evaluating whether genomic features predict the rate of change in anubis allele frequencies over time (Supplementary Methods 14.2).

Briefly, we model the trajectory of anubis ancestry over time for each 100 kb window of the genome as individuals enter or leave the study population, retaining the chance in ancestry per year. We then ask whether genomic features (e.g. recombination, B values, or Fst) predict variation in the change over time across the genome. 

```console
./1anubis_ancestry_change_over_time.R
## For 100 kb windows of the genome, this script requires ancestry information per individual, B values, and recombination rates which are generated in Section 04 and available here as `100kb_ancestry_and_features.RData`. It also requires demographic data (see below) and estimates of Fst per 100 kb window (see below). 
## 
```

#### Data used in these analyses

* **100kb_ancestry_and_features.RData**: contains 100 kb windows of the genome with ancestry calls per individual (`window_by_individual`) and genomic features generated in Section 04 (`new_windows` with columns for B values, `B`, and recombination rate, `rcr`). This file also includes the maximum recombination rate to consider (`max_rcr`), which is 100x the median recombination rate for 100 kb windows of the genome. 
* **amboseli.demographic.info.txt**: contains demographic data for Amboseli individuals with confirmed ids, including the first and last years they were present in the population. Note that other demographic data used in other scripts are also included in this data frame but which we won't use here (their sex, entry type (B = born into a study group, O = first observed in a new study group, or I = immigrating into a study group), genome_wide_anubis_ancestry)
* **FST**: to get these estimates, use the script 4bFST_SNPRCref.sh in "VilgalysFogel_Amboseli_admixture/02_Local ancestry calling/02.3.pedigree_trios/" but substituting 100000 for window size --fst-window-size in place of 35000.
