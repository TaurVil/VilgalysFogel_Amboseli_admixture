## Landscape of introgression

In other parts of this Github (Sections 5 & 6), we show that genomic features are related to the extent of introgressed anubis ancestry in the Amboseli baboons. Here, we take this a step further and ask how much of the total variation in introgressed anubis ancestry can be predicted from genomic features (Section 07.1) using a machine learning approach (elastic net regression). We also ask whether these same genomic features can be used to predict the change in introgressed anubis ancestry over time (Section 07.2)

#### 07.1: using glmnet to predict anubis introgression across the genome

This sub folder contains the data, scripts, and results of using machine learning approaches (elastic net regression using glmnet) to test how much of the total variation in introgressed anubis ancestry can be predicted using genomic features. We note that several features in our dataset are correlated with one another; thus the predictors used produce the best predictive power but it is difficult to interpret which correlated variable directly drives the patterns we observe. 

```console 
## step 1: integrate genomic features matrix
## we expand upon the genomic features available in Section 5 and VilgalysFogel_main_data_file.RData, to include other possible predictors of anubis ancestry such as gene and enhancer content. 

./run.01.get_genomic_features.R
## this script produces ancestry_and_features.RData. It calls the main data file, files from the resources folder, and publicly available baboon annotations from NCBI. 

## step 2: run glmnet
## in R, load the previous matrix and fit elastic net regression using the R package glmnet (Friedman, Hastie, & Tibshirani, 2010). 
## we'll fit 200 iterations, each of which is fit to a random 75% of the genome. We'll repeat this for permuted ancestry values and for recent and historic ancestry. 

./run.02.run_glmnet.R
## this script produces a summary file called glmnet_model_results.RData. It also saves a data file (glmnet_windows.txt) that is presented as Supplementary Table 5. 

## steps 3 and 4: analyze resulting models
./run.03.glmnet_results.R
## a summary of the glmnet results produced by 02.run_glmnet.R. This file includes statistics included in the text and outputs the data table for Fig 4C (glmnet_results.txt). 

./figure4D.R
## produces Fig 4D
```

#### 07.2: change in anubis allele frequencies over time

This directory contains an R script and data files for evaluating whether genomic features predict the rate of change in anubis allele frequencies over time (Supplementary Methods 14.2). Briefly, we first model the trajectory of anubis ancestry across several decades for 100 kb windows along the genome as individuals enter or leave the population. Using the parameter estimate for the year effect as the estimated annual change in anubis ancestry, we then ask whether genomic features (e.g., recombination rate, B values, FST) predict variation in the rate of ancestry change across the genome. 

```console
./run.01.anubis_ancestry_change_over_time.R
## For 100 kb windows of the genome, this script requires ancestry information per individual, B values, and recombination rates which are generated in Sections 04 and 05, available here as `100kb_ancestry_and_features.RData`. It also requires demographic data (see below) and estimates of FST per 100 kb window (see below). 
```

This analysis uses 3 data files:
* **100kb_ancestry_and_features.RData**: contains 100 kb windows of the genome with ancestry calls per individual (`window_by_individual`) and genomic features generated in Section 04 (`new_windows` with columns for B values, `B`, and recombination rate, `rcr`). This file also includes the maximum recombination rate to consider (`max_rcr`), which is 100x the median recombination rate for 100 kb windows of the genome. 
* **amboseli.demographic.info.txt**: contains demographic data for Amboseli individuals with confirmed ids, including the first and last years they were present in the population. Note that other demographic data used in other scripts are also included in this data frame but which we won't use here (their sex, entry type (B = born into a study group, O = first observed in a new study group, or I = immigrating into a study group), genome_wide_anubis_ancestry)
* **FST**: to get these estimates, use the script `4bFST_SNPRCref.sh` in "VilgalysFogel_Amboseli_admixture/02_local-ancestry-calling/02.3.pedigree_trios/" but substituting 100000 for window size --fst-window-size in place of 35000.

