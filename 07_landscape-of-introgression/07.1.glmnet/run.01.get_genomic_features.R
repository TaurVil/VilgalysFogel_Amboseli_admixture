## Prep data
############################

############################
## Load in data
############################
load(paste("./VilgalysFogel_main_data_file.250kb_windows.RData", sep=""))
distance='250kb'

## clean up some unnecessary pieces
rm(d_name, distance, masking, max_rcr, n, unfiltered_for_recombination_rate, ids) 

############################
## get GC and CpG content
############################
gc <- read.delim("Resources/gc_content.panubis1.250kb.txt", header=T)
colnames(gc) <- c("chr", "start", "end", "perc_AT", "perc_GC", "num_A", "num_C", "num_G", "num_T", "num_N", "num_other", "seq_length" )
gc$pos <- paste(gc$chr, gc$start, gc$end, sep="_")

cpg <- read.delim("Resources/cpg_islands.panubis1.emboss.bed", header=F)

## integrate GC and CpG content
colnames(cpg) <- c('chr', 'start', 'end')
gc$islands <- gc$islands_bp <- NULL; cpg$l <- cpg$end - cpg$start
cpg$chr <- as.character(cpg$chr)
for (i in 1:nrow(gc)) {
  subset(cpg, cpg$chr == gc$chr[i] & cpg$start < gc$end[i] & cpg$end > gc$start[i]) -> tmp 
  tmp$start[tmp$start < gc$start[i]] <- gc$start[i]
  tmp$end[tmp$end > gc$end[i] ] <- gc$end[i]
  nrow(tmp) -> gc$islands[i] 
  sum(tmp$l)/distance -> gc$islands_bp[i]
}; rm(i, tmp, cpg)

## integrate with rest of the data frame
as.data.frame(to_analyze) -> data; rm(to_analyze)
data$pos <- paste(data$chr, as.character(data$start+1), sep="_")
gc$pos <- paste(gc$chr, as.character(gc$start +1), sep="_")
subset(data, data$pos %in% gc$pos) -> data
sum(gc$pos == data$pos)/nrow(data); nrow(gc) == nrow(data)
cbind(data, gc$perc_GC, gc$islands, gc$islands_bp) -> data ; rm(gc)

############################
## get gene content
############################
## already have a set of protein coding genes from within section 06
load("Windows_Genes_empty.RData")
## have the results of differentially expression analyses, select for significant effects
mash_res <- read.delim("mash_results.txt")
sig <- mash_res[mash_res$lfsr_eLife < 0.1 | mash_res$lfsr_TC < 0.1,]
empty_windows -> genes; rm(empty_windows)

## convert into 250 kb windows and add to emerging data frame
data$gene_DE <- data$gene_DE_length <- data$gene_exp_in_blood <- data$gene_exp_in_blood_length <- data$gene_number <- data$gene_length <- NA
for (i in 1:nrow(data)) {
  tmp <- subset(genes, genes$chr == data$chr[i] & genes$initial_start < data$end[i] & genes$initial_end > data$start[i]) 
  length(unique(tmp$gene_id)) -> data$gene_number[i]
  if (nrow(tmp) >= 1) {
    possible <- NULL; for (k in 1:nrow(tmp)) {c(possible, tmp$initial_start[k]:tmp$initial_end[k]) -> possible}
    data$gene_length[i] <- sum(data$start[i]:data$end[i] %in% possible); rm(possible)
  }
  
  subset(tmp, tmp$gene_id %in% mash_res$gene) -> tmp
  length(unique(tmp$gene_id)) -> data$gene_exp_in_blood[i]
  if (nrow(tmp) >= 1) {
    possible <- NULL; for (k in 1:nrow(tmp)) {c(possible, tmp$initial_start[k]:tmp$initial_end[k]) -> possible}
    data$gene_exp_in_blood_length[i] <- sum(data$start[i]:data$end[i] %in% possible); rm(possible)
  }
  
  subset(tmp, tmp$gene_id %in% sig$gene) -> tmp 
  length(unique(tmp$gene_id)) -> data$gene_DE[i]
  if (nrow(tmp) >= 1) {
    possible <- NULL; for (k in 1:nrow(tmp)) {c(possible, tmp$initial_start[k]:tmp$initial_end[k]) -> possible}
    data$gene_DE_length[i] <- sum(data$start[i]:data$end[i] %in% possible); rm(possible)
  }
}; rm(i, tmp)
data$yes_no_blood <- 0; data$yes_no_blood[data$gene_exp_in_blood > 0] <- 1
data$perc_DE <- data$gene_DE/data$gene_exp_in_blood

############################
## get enhancer content
############################
enh <- read.delim("~/my_genomes/panubis1/H3K4me1_peaks.bed", header=F) ; colnames(enh) <- c("chr", "start", "end", "enh_name")
data$enh <- NA; for (i in 1:nrow(data)) {
  tmp <- subset(enh, enh$chr == data$chr[i] & enh$start < data$end[i] & enh$end > data$start[i]) 
  tmp$start[tmp$start < data$start[i]] <- data$start[i]
  tmp$end[tmp$end > data$end[i]] <- data$end[i]
  
  tmp$length <- tmp$end - tmp$start
  data$enh[i] <- sum(tmp$length)
  
}; rm(i, tmp, enh)

############################
## add transformed versions of predictors
############################

data$logrcr <- log(data$rcr)
data$rcr_top_quarter <- 0
data$rcr_top_quarter[data$rcr > quantile(data$rcr, 0.75)] <- 1
data$rcr_top_decile <- 0
data$rcr_top_decile[data$rcr > quantile(data$rcr, 0.9)] <- 1

data$B <- (data$B_exons/1000*data$B_regulatory/1000)*1000
data$rankrcr <- rank(data$rcr)
data$perc_DE[is.na(data$perc_DE)] <- 0
data$gene_exp_in_blood_length[is.na(data$gene_exp_in_blood_length)] <- 0
data$gene_DE_length[is.na(data$gene_DE_length)] <- 0
data$presence_DE <- 0; data$presence_DE[data$gene_DE >= 1] <- 1

############################
## cleanup & save RData object
############################
rm(sig, mash_res, distance, genes, k, ped)

data -> including_windows
colnames(data)[24:26] <- c('perc_GC', 'num_cpg_islands', 'perc_cpg_islands')
data[,-c(1:3, 23)] -> data # remove window labels
data$gene_length[data$gene_number == 0] <- 0

save.image("./ancestry_and_features.RData")
