
## step 1: integrate genomic features matrix
## we expand upon the genomic features available in Section 5 and VilgalysFogel_main_data_file.RData, to include other possible predictors of anubis ancestry such as gene and enhancer content. 

./01.get_genomic_features.R
## this script produces 

## step 2: run glmnet
## in R, load the previous matrix and fit elastic net regression using the R package glmnet (Friedman, Hastie, & Tibshirani, 2010). 

./02.run_glmnet.R
## this script 

## step 3: plot resulting values
run figure4C.R
