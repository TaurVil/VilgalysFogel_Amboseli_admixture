## Selection against regulatory divergence

In this section, we will use previously generated expression data from the same Amboseli baboons sequenced in this study to test whether signatures of selection against admixture are stronger when there is evidence of gene regulatory divergence between baboon species. We find: (i) widespread evidence of regulatory divergence between baboon species, (ii) that the effect of ancestry on gene expression is almost always additive and better predicted by local than global ancestry (suggesting cis-regulatory divergence rather than trans effects), and (iii) that there is evidence for stronger selection in regions of the genome where local genetic ancestry affects gene expression.  

Additional background on the expression datasets used in this study can be found in Tung et al. 2015 _eLIFE_, Lea et al. 2018 _PNAS_, and Anderson et al. 2021 _bioRxiv_. eLife dataset refers to samples published in Tung et al. 2015; TruCulture refers to samples from both Lea et al. 2018 and Anderson et al. 2021. 

#### Model effect of local ancestry on gene expression

We began with expression data summarized as read counts per gene and individual, provided here as `./DATA/expression_data.txt` with accompanying metadata `./DATA/expression_metadata.txt`. For each protein coding gene, we calculate the mean anubis ancestry across that gene for each individual.  

```console
./run.01.assign_gene_ancestry.R
# Calls a file of chromosome lengths (downloaded from NCBI), the Panubis1 gtf file (downloaded from NCBI), Amboseli ancestry calls (Zenodo), and recombination rates per chromosome (03_Resources)
# Produces "genes_ancestry.RData", which contains anubis ancestry per gene and individual, mean anubis ancestry per gene across the population, and mean recombination rate per gene
```

Next we process the expression data separately for each dataset. We make the results of this step available as `expression_data_for_models.DATASET.RData`, although the analyses can be reproduced from the count and metadata provided in `/DATA`. 

```console
./run.02.process_expression.R ## run with dataset="eLife" and dataset="TruCulture"
```

Linear models for the effect of local and global genetic ancestry on gene expression. We find an effect of local but not global ancestry on variation in gene expression. We explore the local ancestry effect in more detail below.

```console
./run.03.model_expression.R ## run with dataset="eLife" and dataset="TruCulture"
# Includes 50 permutations of local and global ancestry
# Calls "expression_data_for_models.DATASET.RData" which is returned by run.02.process_expression.R and included in DATA. 
# Outputs a data file "expression_results.DATASET.RData" which is not included here due to file size constraints but can be generated using this script or is available upon request. 
./run.03b.global_ancestry_results.R ## integrate expression data for global ancestry effects, which finds no significant associations. Local ancestry effects, where there is an enrichment, will be integrated below in run.05.mash_linear_models.R.
```

Using piece-wise regression, we identify little evidence for non-additive effects, analyzing all genes where there are at least 10 individuals with each of homozygous yellow, heterozygous, and homozygous anubis ancestry. 

```console 
./run.04.model_expression_nonadditive.R ## run with dataset="eLife" and dataset="TruCulture"
# Calls "expression_data_for_models.DATASET.RData" which is returned by run.02.process_expression.R and included in DATA. 
# Produces "expression_piecewise.DATASET.RData" which is used in the next script
./run.04b.nonadditive_results.R ## summarize model results, using results from run.04.model_expression_nonadditive.R 
```

Because we analyzed two data sets independently, we used multivariate adaptive shrinkage (Urbut et al. 2019 Nature Genetics) to refine effect size estimates.

```console 
## Apply mash to independent analyses of each dataset
./run.05.mash_linear_models.R ## integrate expression data from the two datasets, and refine effect size estimates
# Calls "expression_results.DATASET.RData" from ./run.03.model_expression.R
# Produces statistics for the main text, SI, and "mash_results.txt"
```

#### Selection against genes with greater regulatory divergence

Having observed that local genetic ancestry is associated with expression variation in the Amboseli baboon population, we then asked whether signatures of selection were stronger near genes where ancestry had a larger effect on expression. We find evidence that local introgressed anubis ancestry is reduced in genes where local ancestry has a large effect on gene expression, and that the relationship between recombination rate and introgression is stronger for those same genes. Together, this evidence supports hypotheses that gene regulatory divergence is selected against following admixture. 

```console
./run.06.selection_against_introgressed_regulatory_variation.R
# Calls "mash_results.txt" and "Genes_ancestry.RData" from earlier in this directory. 
# Produces tables with the reduction of introgressed anubis ancestry in genes with ancestry-associated expression ("reduced_ancestry_in_DE_genes.txt") and Spearman's rho between ancestry and recombination rate for the genes with the highest vs. the lowest effects of ancestry on expression ("bootstrap_rho_results.txt").

## Produce figure 4A and 4B
./figure4A_4B.R
# Calls expression data for an example gene, the distribution of p-values, and bootstrap_rho_results.txt to produce the first two panels of Fig 4
```

_Several intermediate files are not included here due to file size constraints, but are available upon request._