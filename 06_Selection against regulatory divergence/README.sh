
## For background on the expression datasets used in this study, see Tung et al. 2015, Lea et al. 2019, and Anderson et al. 2021

## We began with expression data summarized as read counts per gene and individual, provided here as "expression_data.txt" with accompanying metadata "expression_metadata.txt"

## For each protein coding gene, we calculate the ancestry of each gene, we calculate the mean anubis ancestry per individual in the population. 
./01.assign_gene_ancestry.R
# Calls a file of chromosome lengths (downloaded from NCBI), the Panubis1 gtf file (downloaded from NCBI), Amboseli ancestry calls (Zenodo), and recombination rates per chromosome (03_Resources)
# Produces "genes_ancestry.RData", which also contains the mean anubis ancestry and the mean recombination rate per gene

## Script used to process the expression data per dataset. We make the results of this step available as "expression_data_for_models.DATASET.RData", although the analyses can be reproduced from the count and meta data provided in the sub folder within this directory. 
./02.process_expression.R ## run with dataset="eLife" and dataset="TruCulture"

## Run linear models
./03.model_expression.R ## run with dataset="eLife" and dataset="TruCulture"
# Includes 50 permutations of local and global ancestry
# Calls "expression_data_for_models.DATASET.RData" which is returned in step 2 and included in this directory. 
# Outputs a data file "expression_results.DATASET.RData" which is not included here due to file-size considerations but can be generated using this script or is available upon request. 
./03b.global_ancestry_results.R ## integrate expression data for global ancestry effects, which finds no significant associations. Local ancestry effects, where there is an enrichment, will be integrated in step 5.

## Test for non-additive effects
## Note we only analyze data for which there is at least 10 indivuals with homozygous and heterozygous ancestry
./04.model_expression_nonadditive.R ## run with dataset="eLife" and dataset="TruCulture"
# Calls "expression_data_for_models.DATASET.RData" which is returned in step 2 and included in this directory. 
# Produces "expression_piecewise.DATASET.RData" which is used in 04b, but not included here due to file size considerations although it can be generated or is available upon request.
./04b.nonadditive_results.R ## summarize model results, using results from 04 

## Apply mash to independent analyses of each dataset
./05.mash_linear_models.R ## integrate expression data from the two datasets, and refine effect size estimates
# Calls "expression_results.DATASET.RData" from step 3
# Produces statistics for the main text, SI, and "mash_results.txt"

## Selection against genes with greater regulatory divergence
./06.selection_against_introgressed_regulatory_variation.R
# Calls "mash_results.txt" and "genes_ancestry.RData" from earlier in this directory. 
# Produces tables with the reduction of anubis ancestry in genes with ancestry-associated expression ("reduced_ancestry_in_DE_genes.txt") and Spearman's rho between ancestry and recombination for the genes with the highest vs the lowest effects of ancestry on expression ("bootstrap_rho_results.txt").

