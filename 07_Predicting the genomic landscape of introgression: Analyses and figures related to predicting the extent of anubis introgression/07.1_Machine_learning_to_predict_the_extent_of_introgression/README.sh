
## step 1: integrate genomic features matrix
## we expand upon the genomic features available in Section 5 and VilgalysFogel_main_data_file.RData, to include other possible predictors of anubis ancestry such as gene and enhancer content. 

./01.get_genomic_features.R
## this script produces ancestry_and_features.RData

## step 2: run glmnet
## in R, load the previous matrix and fit elastic net regression using the R package glmnet (Friedman, Hastie, & Tibshirani, 2010). 
## we'll fit 200 iterations, each of which is fit to a random 75% of the genome. We'll repeat this for permuted ancestry values and for recent and historic ancestry. 

./02.run_glmnet.R
## this script produces a summary file called glmnet_model_results.RData. It also saves a data file (glmnet_windows.txt) that is presented as Supplementary Table 5. 

## step 3: analyze values

./03.glmnet_results.R
## a summary of the glmnet results produced by 02.run_glmnet.R. This file includes statistics in the text and outputs the data table for Fig. 4C (glmnet_results.txt). 

figure4C.R
## produces fig 4C
