## Get windows for 25kb tracts based on the chromosome lengths (in resources)

## We'll add features, but keep calling it `empty_windows` until the end when it's renamed `windows_25kb`. 
## We'll have a separate object for ancestry: `ancestry_windows`

distance <- 2.5e4
library(data.table); library(parallel)

lengths <- read.delim("../03_baboon-genomic-resources/Resources/Panubis1_chromlengths.txt", header=F)

empty_windows <- NULL 
for (i in 1:20) {
  max <- lengths$V2[i]
  rows <- ceiling(max/distance); print(i); print(rows)
  tmp <- as.data.frame(matrix(ncol=3, nrow=rows))
  colnames(tmp) <- c("chr", "start", "end")
  tmp$chr <- paste("chr",i,sep="")
  tmp$start <- (0:(rows-1))*distance
  tmp$end <- (1:rows)*distance
  tmp <- subset(tmp, tmp$start >= 5e4 & tmp$end <= (max-5e4))
  rbind(empty_windows,tmp) -> empty_windows; rm(tmp, max, rows)
}; rm(i)
rm(lengths)

## add in B values
all2 <- NULL; for (k in 1:20) {
  tmp <- subset(empty_windows, empty_windows$chr == paste("chr",k,sep=""))
  name=paste("../03_baboon-genomic-resources/Resources/B_values_coding/chr",k,".bkgd",sep="")
  b <- fread(name)
  tmp$B_exons <- NA
  if  (nrow(b) > 0) {
    colnames(b) <- c("B", "length"); b$end <- cumsum(b$length); b$start <- c(0,b$end[-nrow(b)])
    for (i in 1:nrow(tmp)) {
      subset(b, b$start <= tmp$end[i] & b$end > tmp$start[i]) -> tmp2
      tmp2$start[tmp2$start < tmp$start[i]] <- tmp$start[i]
      tmp2$end[tmp2$end > tmp$end[i]] <- tmp$end[i]
      tmp2$l <- tmp2$end - tmp2$start
      tmp$B_exons[i] <- sum(tmp2$B*tmp2$l)/sum(tmp2$l); rm(tmp2)
    }; rm(i)
  }; rm(b, name)
  name=paste("../03_baboon-genomic-resources/Resources/B_values_noncoding_near_genes/chr",k,".bkgd",sep="")
  b <- fread(name)
  tmp$B_regulatory <- NA
  if  (nrow(b) > 0) {
    colnames(b) <- c("B", "length"); b$end <- cumsum(b$length); b$start <- c(0,b$end[-nrow(b)])
    for (i in 1:nrow(tmp)) {
      subset(b, b$start <= tmp$end[i] & b$end > tmp$start[i]) -> tmp2
      tmp2$start[tmp2$start < tmp$start[i]] <- tmp$start[i]
      tmp2$end[tmp2$end > tmp$end[i]] <- tmp$end[i]
      tmp2$l <- tmp2$end - tmp2$start
      tmp$B_regulatory[i] <- sum(tmp2$B*tmp2$l)/sum(tmp2$l); rm(tmp2)
    }; rm(i)
  }; rm(b, name)
  all2 <- rbind(all2, tmp); print(k) ; rm(tmp)
}; rm(k)
all2 -> empty_windows; empty_windows$B <- (empty_windows$B_exons/1000*empty_windows$B_regulatory/1000)*1000

## add in recombination rates
all2<- NULL
for (k in 1:20) { #For each chromosome
  # Recombination rates for 24 SNPRC anubis baboons
  name=paste("../03_baboon-genomic-resources/Resources/Recombination_Rates/anubisSNPRC.",k,".txt",sep="")
  fread(name) -> RCR; c(colnames(RCR)[-1],"fill") -> colnames(RCR); rm(name)
  
  rcr_chrom <- subset(empty_windows, empty_windows$chr == paste("chr",k,sep=""))
  rcr_chrom$recombination <- NA
  
  ## calculate mean recombination rate per window
  for (i in 1:nrow(rcr_chrom)) {
    subset(RCR, RCR$left_snp <= rcr_chrom$end[i] & RCR$right_snp > rcr_chrom$start[i]) -> tmp2
    tmp2$left_snp[tmp2$left_snp < rcr_chrom$start[i]] <- rcr_chrom$start[i]
    tmp2$right_snp[tmp2$right_snp > rcr_chrom$end[i]] <- rcr_chrom$end[i]
    tmp2$l <- tmp2$right_snp - tmp2$left_snp
    rcr_chrom$recombination[i] <- sum(tmp2$mean*tmp2$l)/sum(tmp2$l); rm(tmp2)
  }; rm(RCR,i)
  print(c(nrow(rcr_chrom),"chrom", k))
  
  rbind(all2,rcr_chrom) -> all2; rm(rcr_chrom)
}; rm(k)
all2 -> empty_windows

## add in number of variable sites, and number of highly differentiated sites

### Read in allele frequencies and calculate Fst
load("./allele_frequencies.reference.masked.RData")
## get overall yellow and anubis calls rather than split by population
calls[,1:8] -> calls; gc()
calls <- subset(calls, calls$n_all_anubis >= 10 & calls$n_all_yellow >= 10); gc()
paste("chr", calls$chrom, sep="") -> calls$chrom
## remove any sites fixed in both yellow and anubis baboons
calls <- subset(calls, (calls$all_anubis > 0 | calls$all_yellow > 0) & (calls$all_anubis < 1 | calls$all_yellow < 1) )

## calculate fst for each variant
calls$mean_alt_af <- rowMeans(cbind(calls$all_anubis, calls$all_yellow))
ht <- 2*calls$mean_alt_af*(1-calls$mean_alt_af); ha <- 2*calls$all_anubis*(1-calls$all_anubis); hy <- 2*calls$all_yellow*(1-calls$all_yellow)
hsub <- (ha + hy)/2; calls$fst <- (ht-hsub)/ht; rm(ht, hsub, ha, hy); gc()

# hist(calls$fst, xlab="Fst", freq=F, breaks=100)
# hist(calls$fst, xlab="Fst", freq=F)

all2 <- NULL; for (k in 1:20) {
  w2 <- subset(empty_windows, empty_windows$chr == paste("chr",k,sep=""))
  c2 <- subset(calls, calls$chrom == paste("chr", k, sep=""))
  
  for (i in 1:nrow(w2)) {
    tmp2 <- subset(c2, c2$snp <= w2$end[i] & c2$snp > w2$start[i])
    w2$snps[i] <- nrow(tmp2)
    w2$fst50[i] <- sum(tmp2$fst >= 0.5)
    w2$fst55[i] <- sum(tmp2$fst >= 0.55)
    w2$fst60[i] <- sum(tmp2$fst >= 0.60)
    w2$fst65[i] <- sum(tmp2$fst >= 0.65)
    w2$fst70[i] <- sum(tmp2$fst >= 0.70)
    w2$fst75[i] <- sum(tmp2$fst >= 0.75)
    w2$fst80[i] <- sum(tmp2$fst >= 0.80)
    w2$fst85[i] <- sum(tmp2$fst >= 0.85)
    w2$fst90[i] <- sum(tmp2$fst >= 0.90)
    w2$fst95[i] <- sum(tmp2$fst >= 0.95)
    w2$fixed[i] <- sum(tmp2$fst >= 1)
    rm(tmp2)
  }
  all2 <- rbind(all2, w2); print(k); rm(w2,c2,i)
}; rm(k)
all2 -> empty_windows


## Create table of mean ancestry per individual per window
calc_mean_ancestry_per_window <- function(window) {
  subset(t2, t2$start <= walker$end[window] & t2$end > walker$start[window]) -> tmp2
  tmp2$start[tmp2$start < walker$start[window]] <- walker$start[window]
  tmp2$end[tmp2$end > walker$end[window]] <- walker$end[window]
  tmp2$l <- tmp2$end - tmp2$start
  return(sum(tmp2$l*tmp2$state)/sum(tmp2$l))
}

a_data <- fread("./amboseli_LCLAE_tracts.txt"); colnames(a_data) <- c("chr", "start", "end", "state", "length", "u_prev", "u_after", "id")
n <- length(unique(a_data$id))
## remove ancestry tracts less than 1kb long
a_data$length <- a_data$end - a_data$start; a_data <- subset(a_data, a_data$length >= 1000)

ancestry_per_individual <- NULL; 
for (i in 1:20) {
  tmp <- subset(a_data, a_data$chr == paste('chr',i,sep=""))
  win <- subset(empty_windows, empty_windows$chr == paste('chr',i,sep=""))
  
  max_i <- nrow(win)
  walker <- as.data.frame(matrix(ncol=sum(3,n), nrow=max_i) )
  win[,1:3] -> walker[,1:3]; colnames(win)[1:3] -> colnames(walker)[1:3]; rm(win)
  colnames(walker)[4:(n+3)] <- unique(a_data$id)
  
  print(paste("chrom", i)); print(max_i)
  
  for (id in 1:n) {
    nom <- unique(a_data$id)[id]
    t2 <- subset(tmp, tmp$id == nom)
    walker[,id+3] <- do.call("c",mclapply(1:max_i, FUN = calc_mean_ancestry_per_window))
    if (ceiling(id/45) == id/45) {print(paste(nom, "done"))}
  }; rm(id, tmp, win, max_i, t2, nom)
  rbind(ancestry_per_individual,walker) -> ancestry_per_individual; rm(walker)
}; rm(i)

rm(empty_windows, a_data, distance)

rm(calc_mean_ancestry_per_window)

## convert ancestry per individual into the mean ancestry for the population, recent, and historical individuals
## divide the mean by 2, because per individual is the number of alleles
empty_windows -> features

features$mean_ancestry <- rowMeans(ancestry_per_individual)/2
tmp <- ids$table_s1_id[ids$recent_hybridsc == 'historical']
features$historical_ancestry <- rowMeans(ancestry_per_individual[,colnames(ancestry_per_individual) %in% tmp])/2
## for recent ancestry, only include individuals with >40% anubis ancestry on that chromosome
f2 <- NULL; for (chrom in unique(features$chr)) {
  tmp <- rownames(recent)[recent[,chrom] > 0.4]
  tmp_features <- features[features$chr == chrom,]
  tmp_ancestry <- ancestry_per_individual[features$chr == chrom,]
  tmp_features$recent_ancestry <- rowMeans(tmp_ancestry[, colnames(tmp_ancestry) %in% tmp])/2
  rbind(f2, tmp_features) -> f2; rm(tmp_features)
  rm(tmp_ancestry, tmp)
}; rm(chrom)
f2 -> features; rm(f2, ancestry_per_individual)

## Save R object for future steps


save.image("./windows.25kb.RData")
