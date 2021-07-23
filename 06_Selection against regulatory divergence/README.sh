
## For background on the expression datasets used in this study, see Tung et al. 2015, Lea et al. 2019, and Anderson et al. 2021

## We began with expression data summarized as read counts per gene and individual, provided here as "expression_data.txt" with accompanying metadata "expression_metadata.txt"

## For each protein coding gene, we calculate the ancestry of each gene, we calculate the mean anubis ancestry per individual in the population. 
./01.assign_gene_ancestry.R
# Calls a file of chromosome lengths (downloaded from NCBI), the Panubis1 gtf file (downloaded from NCBI), Amboseli ancestry calls (Zenodo), and recombination rates per chromosome (03_Resources)
# Produces genes_ancestry.RData, which also contains the mean anubis ancestry and the nean recombination rate per gene

## We process each dataset to normalize and control for covariates. 
./02.process_expression.R ## run with dataset="eLife" and dataset="TruCulture"

## Run linear models

## Test for non-additive effects

## Apply mash to independent analyses of each dataset

## Selection against genes with greater regulatory divergence
