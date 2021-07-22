
## For background on the expression datasets used in this study, see Tung et al. 2015, Lea et al. 2019, and Anderson et al. 2021

## We began with expression data summarized as read counts per gene and individual, provided here as "expression_data.txt" with accompanying metadata "expression_metadata.txt"

## We process each dataset to normalize and control for covariates. 
./01.process_expression.R ## run with dataset="eLife" and dataset="TruCulture"

## For each protein coding gene, we calculate the ancestry of each gene, we calculate the ancestry 
./02.assign_gene_ancestry.R

## Run linear models

## Test for non-additive effects

## Apply mash to independent analyses of each dataset

## Selection against genes with greater regulatory divergence
