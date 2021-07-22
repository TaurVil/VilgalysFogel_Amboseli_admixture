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
read.delim("~/read_length_per_file.txt", header=F) -> file_read_lengths ## files are either 126 or 101 bp. When in doubt, treat files as 100 bp long. This is only used to calculate TPM and remove lowly expressed genes. 
tpm <- data; for (i in 1:ncol(tpm)) {tmp_length <- file_read_lengths[grepl(pattern = info$ids[i], x=file_read_lengths$V1),]
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

#######################################
## Save output file
#######################################

save.image(paste("./expression_data_", dataset, ".RData", sep=""))
