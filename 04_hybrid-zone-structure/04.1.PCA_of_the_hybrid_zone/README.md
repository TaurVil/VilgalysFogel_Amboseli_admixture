
## Principal Component Analysis

To visualize the genetic variation in our sample, we used principal component analyses to visualize the major axes of variation. We include only high coverage Amboseli genomes (n=9) along with all possible reference panel individuals. We extract these samples from the merged vcf including both Amboseli and reference panel individuals, after filtering (`merged_shared.vcf.gz`, see Section 01_Genotype_Calling), and generate a 012 genotype matrix from which to estimate the covariance between samples. 

```console
## filter for select individuals, and output as 012 genotype matrix
vcftools --gzvcf merged_shared.allchroms.vcf.gz --keep 00_anu.list --keep 00_yel.list --keep 00_hicov_amboseli.list --012 --out for_pca
```

In R, we will then process these data to get the covariance between samples and run PCA. 
```console
library(data.table)
fread("./for_pca.012") -> d; t(d) -> d; d[-1,] -> d; d[d==-1] <- NA

fread("./for_pca.012.indv", head=F) -> names
names$missing <- colSums(is.na(d))
## Add in the population where each animal originated 
names$source <- c(rep("Amboseli",9), rep("SNPRCanubis_founders",24), rep("Mara", 7), rep("WNPRC",6), rep("BGDPanubis",4), rep("SNPRCyellow_founders",7), rep("Mikumi",15))

missing <- rowSums(is.na(d)); hist(missing) # plot amount of missing data per site

## covariance of all individuals, using sites called in almost all individuals (we will allow 1 sample missing each site)
cov(scale(d[missing<=1,], center=T, scale=T), use="pairwise") -> covgeno

## covariance of high coverage individuals, using sites with no missing data
d2 <- d[,-c(34:46,58:60,65)] # exclude the low coverage individuals
m2 <- rowSums(is.na(d2)) # look at the missing data per row
covgeno_highcov <- cov(scale(d2[m2==0,], center=T, scale=T), use="pairwise")  # get the covariance matrix from the scaled and centered genotype matrix only for rows with no missing data (i.e. m2 is 0)

# remove large files of genotype calls, retaining just the covariance matrices and sample data
rm(d, d2, missing, m2)
save.image("for_pca.RData")
```

Finally, we will produce Figures 1B and S4. The data file needed to reproduce these figures is included here (`for_pca.RData`) and the code is in `figure1B_S4.R`
```console
./figure1B_S4.R
```
