################
# Process expression data
################
library(data.table); library(Rgb)

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

## Read in file of protein coding genes and gene lengths (sum of exonic sequence) profiled in the expression data
## this information was extracted from the NCBI genome annotation for Panubis1: GCF_008728515.1_Panubis1.0_genomic.gtf
## gtf <- read.gtf("./GCF_008728515.1_Panubis1.0_genomic.gtf"); genes <- subset(gtf, gtf$feature == "gene" & gtf$gene_biotype == "protein_coding"); exons <- subset(gtf, gtf$feature == "exon" & gtf$gene_id %in% genes$gene_id)
## subset(genes, genes$gene_id %in% row.names(expr)) -> genes; genes[order(genes$gene_id), ] -> genes; rm(gtf); gc()
## 
gene_lengths <- as.data.frame(matrix(ncol=2, nrow=nrow(expression_data))); genes$gene_id -> gene_lengths[,1]
for (i in 1:nrow(gene_lengths)) { subset(exons, exons$gene_id %in% gene_lengths[i,1]) -> tmp;  sum(tmp$end-tmp$start)/length(unique(tmp$transcript_id)) -> gene_lengths[i,2] }; rm(i, tmp)
colnames(gene_lengths)<-c("gene_ID","length"); rm(genes, exons); write.table(gene_lengths, "./gene_lengths.txt", row.names=F, col.names=T, sep="\t", quote=F)

expression_data -> data


## Restrict dataset to protein coding genes
subset(expr, row.names(expr) %in% genes$gene_id) -> expression_data
subset(genes, genes$gene_id %in% row.names(expression_data)) -> genes; 

## Order genes and expression data alphabetically
genes[order(genes$gene_id), ] -> genes
expression_data[order(row.names(expression_data)),] -> expression_data
# sum(genes$gene_id == row.names(expression_data)) # check that the order is the same for genes and the expression data
