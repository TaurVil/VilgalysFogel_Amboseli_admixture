## Selection against regulatory divergence

In this section, we will use previously generated expression data from the same Amboseli baboons sequenced in this study to test whether signatures of selection against admixture are stronger when there is evidence of gene regulatory divergence between baboon species. We find: (i) widespread evidence of regulatory divergence between baboon species, (ii) that the effect of ancestry on gene expression is almost always additive and better predicted by local than global ancestry (suggesting cis-regulatory divergence rather than trans effects), and (iii) that there is evidence for stronger selection in regions of the genome where local genetic ancestry affects gene expression.  

Additional background on the expression datasets used in this study can be found in Tung et al. 2015 _eLIFE_, Lea et al. 2019 _PNAS_, and Anderson et al. 2021 _bioRxiv_. 

#### Model effect of local ancestry on gene expression

We began with expression data summarized as read counts per gene and individual, provided here as "expression_data.txt" with accompanying metadata "expression_metadata.txt"

For each protein coding gene, we calculate the mean anubis ancestry across that gene for each individual.  

```console
./r01.assign_gene_ancestry.R
# Calls a file of chromosome lengths (downloaded from NCBI), the Panubis1 gtf file (downloaded from NCBI), Amboseli ancestry calls (Zenodo), and recombination rates per chromosome (03_Resources)
# Produces "genes_ancestry.RData", which also contains the mean anubis ancestry and the mean recombination rate per gene
```

Next we process the expression data separately for each dataset. We make the results of this step available as `expression_data_for_models.DATASET.RData`, although the analyses can be reproduced from the count and metadata provided in `/DATA`. 

```console
./r02.process_expression.R ## run with dataset="eLife" and dataset="TruCulture"
```

Linear models for the effect of local and global genetic ancestry on gene expression. Local ancestry is explored in more detail later, while global ancestry does not predict variation in gene expression when the model also includes local ancestry. 

```console
./r03.model_expression.R ## run with dataset="eLife" and dataset="TruCulture"
# Includes 50 permutations of local and global ancestry
# Calls "expression_data_for_models.DATASET.RData" which is returned in step 2 and included in this directory. 
# Outputs a data file "expression_results.DATASET.RData" which is not included here due to file size constraints but can be generated using this script or is available upon request. 
./r03b.global_ancestry_results.R ## integrate expression data for global ancestry effects, which finds no significant associations. Local ancestry effects, where there is an enrichment, will be integrated in step 5.
```

Using piece-wise regression, we identify little evidence for non-additive effects, analyzing all genes where there are at least 10 individuals with homozygous yellow, heterozygous, and homozygous anubis ancestry

```console 
./r04.model_expression_nonadditive.R ## run with dataset="eLife" and dataset="TruCulture"
# Calls "expression_data_for_models.DATASET.RData" which is returned in step 2 and included in this directory. 
# Produces "expression_piecewise.DATASET.RData" which is used in 04b, but not included here due to file size constraints although it can be generated using this script or is available upon request.
./r04b.nonadditive_results.R ## summarize model results, using results from 04 
```

Because we analyzed two data sets independently, we used multivariate adaptive shrinkage (Urbut et al. 2019 Nature Genetics) to refine effect size estimates.

```console 
## Apply mash to independent analyses of each dataset
./r05.mash_linear_models.R ## integrate expression data from the two datasets, and refine effect size estimates
# Calls "expression_results.DATASET.RData" from step 3
# Produces statistics for the main text, SI, and "mash_results.txt"
```

#### Selection against genes with greater regulatory divergence

Having observed that local genetic ancestry is associated with expression variation in the Amboseli baboon population, we then asked whether signatures of selection were stronger . We find evidence that local introgressed anubis ancestry is reduced in genes where local ancestry has a large effect on gene expression, and that the relationship between recombination rate and introgression is stronger for those same genes. Together, this evidence supports hypotheses that gene regulatory divergence is selected against following admixture. 

```console
./r06.selection_against_introgressed_regulatory_variation.R
# Calls "mash_results.txt" and "genes_ancestry.RData" from earlier in this directory. 
# Produces tables with the reduction of introgressed anubis ancestry in genes with ancestry-associated expression ("reduced_ancestry_in_DE_genes.txt") and Spearman's rho between ancestry and recombination rate for the genes with the highest vs. the lowest effects of ancestry on expression ("bootstrap_rho_results.txt").

## Produce figure 4A and 4B
./figure4A_4B.R
# Calls expression data for an example gene, the distribution of p-values, and bootstrap_rho_results.txt to produce the first two panels of Fig 4
```
