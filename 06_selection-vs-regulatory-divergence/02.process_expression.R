################
# Process expression data
################
library(data.table); library(Rgb); library(edgeR); library(limma); library(impute)

## Specify dataset
dataset="TruCulture" 

## Read in data
expr <- read.delim("./expression_data.txt", header=T)
## Read in metadata
info <- read.delim("./expression_metadata.txt", header=T) 
## Read in sample metadata from WGS
ids <- read.delim("./amboseli_ids.txt", header=T)

### SELECT DATASET ###
expr[,info$dataset == dataset] -> expr
info[info$dataset == dataset,] -> info

## Remove expression and metadata for an individual which seems to be an incorrect id
expr <- expr[,!(info$studyid == 'AMB_210')] 
info <- info[!(info$studyid == 'AMB_210'),] 

## Restrict to only the individuals whom we have local ancestry and expression data for
expr[,info$studyid %in% ids$table_s1_id] -> expr
info[info$studyid %in% ids$table_s1_id, ] -> info
ids[ids$table_s1_id %in% info$studyid, ] -> ids

## Add mean ancestry into the metadata from the ids table
info$mean_ancestry <- NA; for (i in 1:nrow(info)) {
  subset(ids, ids$table_s1_id == as.character(info$studyid[i])) -> tmp
  if (nrow(tmp) == 1) {info$mean_ancestry[i] <- tmp$genome_wide_anubis_ancestryd }
}; rm(i, tmp)

## Samples are ordered by "Order" in the info column
colnames(expr) <- paste("order",info$Order,sep="")
info$Order2 <- paste("order",info$Order,sep="")

## Read in file of protein coding genes and gene lengths (sum of exonic sequence) profiled in the expression data
## this information was extracted from the NCBI genome annotation for Panubis1: GCF_008728515.1_Panubis1.0_genomic.gtf
## gtf <- read.gtf("./GCF_008728515.1_Panubis1.0_genomic.gtf"); genes <- subset(gtf, gtf$feature == "gene" & gtf$gene_biotype == "protein_coding"); exons <- subset(gtf, gtf$feature == "exon" & gtf$gene_id %in% genes$gene_id)
## subset(genes, genes$gene_id %in% row.names(expr)) -> genes; genes[order(genes$gene_id), ] -> genes; rm(gtf); gc()
## gene_lengths <- as.data.frame(matrix(ncol=2, nrow=nrow(expression_data))); genes$gene_id -> gene_lengths[,1]
## for (i in 1:nrow(gene_lengths)) { subset(exons, exons$gene_id %in% gene_lengths[i,1]) -> tmp;  sum(tmp$end-tmp$start)/length(unique(tmp$transcript_id)) -> gene_lengths[i,2] }; rm(i, tmp)
## colnames(gene_lengths)<-c("gene_ID","length"); rm(genes, exons); write.table(gene_lengths, "./gene_lengths.txt", row.names=F, col.names=T, sep="\t", quote=F)
read.delim("./gene_lengths.txt", header=T) -> gene_lengths

## Restrict dataset to protein coding genes
subset(expr, row.names(expr) %in% gene_lengths$gene_ID) -> data ; rm(expr)

## Order genes and expression data alphabetically
gene_lengths[order(gene_lengths$gene_ID), ] -> gene_lengths
data[order(row.names(data)),] -> data
# sum(gene_lengths$gene_ID == row.names(data)) # check that the order is the same for genes and the expression data

#######################################
## Calculate RPKM and TPM
## Filter for genes with TPM >= 2 and measure in >1/2 the dataset
#######################################
info$reads <- colSums(data)

#RPKM
rpkm<-10^9*data/gene_lengths$length
rpkm2<-rpkm
for (i in 1:ncol(data)){rpkm2[,i]<-rpkm[,i]/info$reads[i]}; rm(i)
rpkm2 -> rpkm; rm(rpkm2)

#TPM
read.delim("./read_length_per_file.txt", header=T) -> file_read_lengths ## files are either 126 or 101 bp. When in doubt, treat files as 100 bp long. This is only used to calculate TPM and remove lowly expressed genes. 
tpm <- data; for (i in 1:ncol(tpm)) {tmp_length <- file_read_lengths[grepl(pattern = info$studyid[i], x=file_read_lengths$V1),]
  if (nrow(tmp_length) > 0) {tmp_length <- mean(tmp_length$V2)-1} else {tmp_length <- 100}
  tpm[,i] <- 10^6*data[,i]*tmp_length/gene_lengths$length
}; rm(i, tmp_length, file_read_lengths)

t<-apply(X = data,2,function(x) sum(x*100/gene_lengths$length))
for (i in 1:ncol(data)){
  tpm[,i]<-tpm[,i]/t[i]
}; rm(i, t)

tmp <- data; tmp[tmp==0] <- NA; rowSums(is.na(tmp)) -> missing; rm(tmp)
      
# Filter for sites with mean TPM >= 2 and data in at least half of individuals
subset(data, rowMeans(tpm) >= 2 & missing < ncol(data)/2) -> data
gene_lengths <- subset(gene_lengths, gene_lengths$gene_ID %in% row.names(data)) 
rm(missing, rpkm, tpm)
         
#######################################
## Voom normalization
#######################################
# Create DGE list object
d0 <- DGEList(data)
# Calculate normalization factors
d0 <- calcNormFactors(d0)
## we've already filtered out lowly expressed genes, so we don't need to do so again here
d <- d0; rm(d0)

## create variable for year and dataset
lane <- interaction(info$dataset, info$year)
# plotMDS(d, col = as.numeric(lane))

# Specify the model to be fitted. We do this before using voom since voom uses variances of the model residuals (observed - fitted)
mm <- model.matrix(~1 + info$mean_ancestry + info$age + info$sex)
# The above specifies a model where each coefficient corresponds to a group mean
# Voom
y <- voom(d, mm, plot = F); rm(d, data)
# plotMDS(y$E, col = as.numeric(lane))

#######################################
## Use limma to regress out effects of Batch/year & flow data (for TC)
#######################################
info$yr <- paste("yr", info$year,sep="")
      
# For the truculture dataset, add in the flow cytometry data
# Impute missing values based on the mean value
if (dataset == "TruCulture") {
  metadata2 <- read.delim("./TruCulture_flow_cytometry.txt", header=T)
  t <- merge(x=info, y=metadata2, by = c("studyid", "year"), all.x = T)
  t <- t[order(t$studyid, t$Order),]
  t$yr <- as.factor(t$yr)
  rm(metadata2, info)
  t -> metadata; rm(t) 
  # Impute missing cell composition data
  impute.knn(as.matrix(metadata[,12:16]), rowmax = 1, colmax = 1) -> imp_cellcomp
  # Get PCs of cell composition; colnames(metadata[,34:38])
  pc_comp <- prcomp(imp_cellcomp$data, scale. = T); rm(imp_cellcomp)
  metadata$pc1 <- pc_comp$x[,1]; metadata$pc2 <- pc_comp$x[,2]; metadata$pc3 <- pc_comp$x[,3]; rm(pc_comp)
} else {info -> metadata; rm(info)}
# plotMDS(y$E, col = as.numeric(as.factor(metadata$yr)))
# plotMDS(y$E, col = as.numeric(as.factor(metadata$sex)))
metadata <- metadata[order(metadata$Order),]
         
## use LIMMA to get residuals 
# Construct design matrix to get the residuals
options(na.action='na.pass')
if (dataset == "TruCulture") {design_to_reg_out <- model.matrix(~ 1 + metadata$yr + metadata$pc1 + metadata$pc2)} else {design_to_reg_out <- model.matrix(~ 1 + metadata$yr)}
# Run limma and get residuals
covariate_models <- lmFit(y, design = design_to_reg_out)
residuals <- residuals.MArrayLM(object = covariate_models, y = y$E)
rm(covariate_models, design_to_reg_out)

# plotMDS(residuals, col = as.numeric(as.factor(metadata$yr)))
# plotMDS(residuals, col = as.numeric(as.factor(metadata$sex)))

#######################################
## add in a genetic covariance matrix, called from RNA sequencing data
#######################################
read.delim("./rnaseq_covariance.txt") -> gt; library(data.table)
## identify the individuals and order we want from gt
g <- as.data.frame(gt); row.names(g) <- colnames(g)
for (i in 1:nrow(metadata)) {
  as.data.frame(g[, colnames(g) == metadata$studyid[i]])  -> tmp
  colnames(tmp)[1] <- metadata$Order2[i]
  if (i > 1) {cbind(gt2,tmp) -> gt2}
  if (i == 1) {tmp -> gt2}
}
row.names(gt2) <- row.names(g)
for (i in 1:nrow(metadata)) {
  as.data.frame(gt2[row.names(gt2) == metadata$studyid[i], ])  -> tmp
  rownames(tmp)[1] <- metadata$Order2[i]
  if (i > 1) {rbind(gt3,tmp) -> gt3}
  if (i == 1) {tmp -> gt3}
}
gt3 -> gt_cov2
rm(gt, g, gt2,gt3, tmp, i)
         
#######################################
## add in ancestry per gene, which we generated in step 01
#######################################

load("./Genes_ancestry.RData"); rm(genes, n, features_genes)

## get just the genes we both calculated ancestry and measured expression for
subset(ancestry_genes, ancestry_genes$gene_id %in% row.names(residuals)) -> ancestry; rm(ancestry_genes)
subset(gene_lengths, gene_lengths$gene_ID %in% ancestry$gene_id) -> gene_lengths

## get ancestry matrix for just the target individuals, in the correct order
ancestry -> ancestry_all
for (i in 1:nrow(metadata)) {
  keep <- colnames(ancestry_all) == metadata$studyid[i]
  tmp <- ancestry_all[,keep]
  if (i == 1) {as.data.frame(matrix(ncol=1, data=tmp)) -> ancestry}
  if (i > 1) {cbind(ancestry, tmp) -> ancestry}
  colnames(ancestry)[i] <- metadata$studyid[i]
}; rm(i, tmp, keep)
apply(ancestry, 1, sd) -> sd_ancestry; #hist(sd_ancestry); hist(rowMeans(ancestry, na.rm=T))

## remove genes for which ancestry is largely invariant
row.names(ancestry) <- ancestry_all$gene_id
loc_anc <- subset(ancestry, rowMeans(ancestry, na.rm=T) < 1.8 & rowMeans(ancestry, na.rm=T) > 0.2 & sd_ancestry > 0.4)
loc_anc <- loc_anc[order(row.names(loc_anc)),]
rm(sd_ancestry, ancestry, ancestry_all) # head(loc_anc[,1:7])

## limit to genes where we have ancestry information (removes X genes, which is why the structure changes here to show no sex effects)
subset(gene_lengths, gene_lengths$gene_ID %in% row.names(loc_anc)) -> gene_lengths
subset(residuals, row.names(residuals) %in% row.names(loc_anc)) -> residuals
colnames(loc_anc) <- metadata$Order2

## add design matrix, each individual represented once in the covariance
Z=diag(ncol(loc_anc))
         
#######################################
## cleanup and save image for analyses
#######################################
metadata$sex <- as.factor(metadata$sex)
rm(y, mm, ids)
         
save.image(paste("./expression_data_for_models.", dataset, ".RData", sep=""))

