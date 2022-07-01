## calculate proportion of anubis ancestry per protein coding gene and individual
library(Rgb);library(data.table);library(parallel)

#######################################
## Get properties of genes
#######################################

## read in Panubis1 chromosome lengths (from the fasta index file)
lengths <- read.delim("./Panubis1_chromlengths.txt", header=F)

## read in Panubis1 gtf file
gtf <- read.gtf("./GCF_008728515.1_Panubis1.0_genomic.gtf") # available from NCBI

## identify protein coding genes from the gtf file 
genes <- subset(gtf, gtf$feature == "gene" & gtf$gene_biotype == "protein_coding")
rm(gtf)

## for each chromosome, update the sequence name away from the ncbi strand to 1:20
for (i in 1:20) {genes$seqname[genes$seqname == unique(genes$seqname)[i]] <- i}; rm(i)

# length(unique(genes$gene_id)); length(unique(genes$db_xref)) # get number of genes

# get gene lengths
genes$length <- genes$end - genes$start

#######################################
## create empty windows to populate with information later
#######################################

empty_windows <- genes[,c(1,4:5,7,9)]
colnames(empty_windows)[1] <- 'chr'

empty_windows$initial_start <- empty_windows$start
empty_windows$initial_end <- empty_windows$end
empty_windows$length <- empty_windows$initial_end - empty_windows$initial_start
empty_windows$middle <- rowMeans(cbind(empty_windows$initial_start, empty_windows$initial_end))

#######################################
## for genes less than 50kb long, add sequence centered on the gene to make it 50kb
## this was done because local ancestry and recombination rate are estimated less reliably upon small scales, but does not substantially alter results
#######################################

empty_windows$start[empty_windows$length < 50000] <- empty_windows$middle[genes$length < 50000]-25000
empty_windows$end[empty_windows$length < 50000] <- empty_windows$middle[genes$length < 50000]+25000
empty_windows$start[empty_windows$start < 1] <- 1

# make sure none of the genes are not longer than the chromosomes
empty_windows$start <- apply(cbind(empty_windows$start,1),1,max)
for (i in 1:20) {
  maxlength <- lengths$V2[i] 
  empty_windows$end[empty_windows$chr == i] <- apply(cbind(empty_windows$end[empty_windows$chr == i],maxlength),1,min)
}; rm(i, maxlength)

## adjust chromosome names by adding "chr" prefix
empty_windows$chr <- paste("chr",empty_windows$chr,sep="")

#######################################
## remove genes near the ends of chromosomes (within 50kb) where ancestry can't be estimated accurately
#######################################
new <- NULL; for (i in 1:20) {
  tmp <- subset(empty_windows, empty_windows$start > 50000 & empty_windows$end  < lengths$V2[i] & empty_windows$chr == paste("chr",i,sep="") )
  rbind(new, tmp) -> new; rm(tmp)
}; rm(i)
## reduces dataset from 21091 genes to 19103

new -> empty_windows; rm(new, lengths)
genes <- subset(genes, genes$gene_id %in% empty_windows$gene_id)

## TO REMOVE: called this step "Windows_Genes_empty.RData"
load("~/Baboon/my_walkers/Windows_Genes_empty.RData")

#######################################
## get matrix of individual ancestry per gene : ancestry_genes
#######################################
calc_mean_ancestry_per_window <- function(window) {
  subset(t2, t2$start <= walker$end[window] & t2$end > walker$start[window]) -> tmp2
  tmp2$start[tmp2$start < walker$start[window]] <- walker$start[window]
  tmp2$end[tmp2$end > walker$end[window]] <- walker$end[window]
  tmp2$l <- tmp2$end - tmp2$start
  return(sum(tmp2$l*tmp2$state)/sum(tmp2$l))
}

# load ancestry tracts 
tracts <- fread("./amboseli.tracts.maskedSNPRCref.txt")
# remove tracts <1 kb in length
tracts <- tracts[tracts$length2>1000,]

## get number of samples in the dataset
n <- length(unique(tracts$table_s1_id))

## per chromosome, then per individual, get the mean anubis ancestry for that individual per gene
ancestry_genes <- NULL; 
for (i in 1:20) {
  tmp <- subset(tracts, tracts$chrom == paste('chr',i,sep=""))
  win <- subset(empty_windows, empty_windows$chr == paste('chr',i,sep=""))
  
  max_i <- nrow(win)
  walker <- as.data.frame(matrix(ncol=sum(5,n), nrow=max_i) )
  win[,1:5] -> walker[,1:5]; colnames(win)[1:5] -> colnames(walker)[1:5]; rm(win)
  colnames(walker)[6:(n+5)] <- unique(tracts$table_s1_id)
  
  print(paste("chr", i)); print(max_i) # report chromosome and number of genes on that chromosome
  
  ## step per individual
  for (id in 1:n) {
    nom <- unique(tracts$table_s1_id)[id]
    t2 <- subset(tmp, tmp$table_s1_id == nom)
    walker[,id+5] <- do.call("c",mclapply(1:max_i, FUN = calc_mean_ancestry_per_window))
    if (ceiling(id/45) == id/45) {print(paste(nom, "done"))} #  report progress
  }; rm(id, tmp, max_i, t2, nom)
  
  rbind(ancestry_genes,walker) -> ancestry_genes; rm(walker)
}; rm(i)

rm(tracts, calc_mean_ancestry_per_window)

#######################################
## get recombination per gene
#######################################
rcr<- NULL
for (k in 1:20) { #For each chromosome
  # n24 SW anubis recombination rates. These values can be found on GitHub in the 03_Resources folder under Recombination
  name=paste("./anubisSNPRC.",k,".txt",sep="")
  fread(name) -> RCR; c(colnames(RCR)[-1],"fill") -> colnames(RCR); rm(name)
  
  rcr_chrom <- subset(empty_windows, empty_windows$chr == paste("chr",k,sep=""))
  rcr_chrom$n24_anubis <- NA
  
  for (i in 1:nrow(rcr_chrom)) {
    
    subset(RCR, RCR$left_snp <= rcr_chrom$end[i] & RCR$right_snp > rcr_chrom$start[i]) -> tmp2
    tmp2$left_snp[tmp2$left_snp < rcr_chrom$start[i]] <- rcr_chrom$start[i]
    tmp2$right_snp[tmp2$right_snp > rcr_chrom$end[i]] <- rcr_chrom$end[i]
    tmp2$l <- tmp2$right_snp - tmp2$left_snp
    rcr_chrom$n24_anubis[i] <- sum(tmp2$mean*tmp2$l)/sum(tmp2$l); rm(tmp2)
  }; rm(RCR,i)
  
  print(c(nrow(rcr_chrom),"chrom", k))
  
  rbind(rcr,rcr_chrom) -> rcr; rm(rcr_chrom)
}; rm(k)
empty_windows$recombination_rate <- rcr$n24_anubis; rm(rcr)

#######################################
## get mean population-wide anubis ancestry per gene, as well as recent anubis ancestry
#######################################
load("./tmp_genes.RData")

## briefly load in the overall data frame to get recent/historical individuals
load("./VilgalysFogel_main_data_file.250kb_windows.RData"); rm(to_analyze, unfiltered_for_recombination_rate, d_name, distance, masking, max_rcr)

## mean ancestry
empty_windows$mean_ancestry <- rowMeans(ancestry_genes[,6:(5+n)], na.rm=T)/2
## historical hybrids
t_hist_names <- ids$table_s1_id[ids$recent_hybridsc == 'historical']
empty_windows$historical_ancestry <- rowMeans(ancestry_genes[,colnames(ancestry_genes) %in% t_hist_names], na.rm=T)/2
rm(t_hist_names)
## recent hybrids, including chromosomes with > 40% anubis ancestry
fun_get_recent <- function(site, input_ancestry) {
  t_recent_names <- rownames(recent)[recent[,empty_windows$chr[site]] > 0.4]
  t_recent <- input_ancestry[site,colnames(input_ancestry) %in% t_recent_names]
  return(mean(t(t_recent), na.rm=T)/2)
}

system.time(r <- mclapply(1:nrow(empty_windows), FUN=fun_get_recent, input_ancestry=ancestry_genes))
r <- do.call("rbind",r); empty_windows$recent_ancestry <- r[,1]
empty_windows -> features_genes; rm(empty_windows, fun_get_recent, recent, r)

#######################################
## save file to be called next
#######################################
save.image("./Genes_ancestry.RData")
