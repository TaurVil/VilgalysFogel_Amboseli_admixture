library(data.table)

# read in chromosome file 
read.delim("~/genomes/panubis1/chromInfo.txt", header=F) -> chroms; colnames(chroms)[1:2] <- c("name", "length")
# read in file of test and source individuals 
read.delim("./DATA/00_yellow_sources.list", header=F) -> yel_source; read.delim("./DATA/00_yel.list", header=F) -> yel
read.delim("./DATA/00_anubis_sources.list", header=F) -> anu_source; read.delim("./DATA/00_anu.list", header=F) -> anu

# read in all yellow tracts
IBD_yellow <- NULL; for (i in 1:nrow(yel_source)) { tmp <- NULL; for (chrom in 1:20) {
    name=paste("./IBDmix_by_chrom/yellow.relative_to_",yel_source[i,1],".",chrom,".txt",sep="")
    fread(name) -> data; data$length <- data$end - data$start; data$source <- yel_source[i,1]; data <- subset(data, data$length >= 1000 & data$slod >= 4)
    rbind(tmp,data) -> tmp }; rbind(IBD_yellow,tmp) -> IBD_yellow; rm(tmp); print(i) }; rm(i, chrom, data, name)
IBD_yellow <- IBD_yellow[order(IBD_yellow$ID,IBD_yellow$source,IBD_yellow$chrom,IBD_yellow$start),]

# read in all anubis tracts
IBD_anubis <- NULL; for (i in 1:nrow(anu_source)) { tmp <- NULL; for (chrom in 1:20) {
    name=paste("./IBDmix_by_chrom/anubis.relative_to_",anu_source[i,1],".",chrom,".txt",sep="")
    fread(name) -> data; data$length <- data$end - data$start; data$source <- anu_source[i,1]; data <- subset(data, data$length >= 1000 & data$slod >= 4)
    rbind(tmp,data) -> tmp }; rbind(IBD_anubis,tmp) -> IBD_anubis; rm(tmp); print(i) }; rm(i, chrom, data, name)
IBD_anubis <- IBD_anubis[order(IBD_anubis$ID,IBD_anubis$source,IBD_anubis$chrom,IBD_anubis$start),]

# convert files to a per-segment readout with how many times each segment is called as IBD in possible source individuals
# for yellow baboons, only use SNPRC anubis baboons and the unadmixed individual from Aberdares as source individuals
# for anubis baboons, only use Mikumi yellow baboons as the source individuals
library(IRanges); library(data.table)
res_collapsed <- NULL; for (i in 1:nrow(yel)) {
  tmp <- IBD_yellow[IBD_yellow$ID == yel$name[i] & IBD_yellow$source %in% yel_source$name[c(1:24,26:28)],]; tmp2 <- NULL
  for (chr in 1:20) {
    t3 <- tmp[tmp$chrom == chr,]; t3 <- t3[order(t3$start, t3$end),]
    cov <- coverage(IRanges(t3$start, t3$end)); as.data.table(cbind(cov@values, cov@lengths)) -> mc
    mc$end <- cumsum(mc$V2); mc$start <- c(1,mc$end[-length(mc$end)]+1); colnames(mc) <- c("number", "length", "end", "start"); mc$chrom <- chr
    tmp2 <- rbind(tmp2[,c(5,4,3,1,2)], mc); rm(mc,cov,t3)
  }; rm(chr)
  tmp2$name <- yel$name[i]; res_collapsed <- rbind(res_collapsed, tmp2); rm(tmp2)
}; rm(i, tmp, IBD_yellow); yel_collapsed <- res_collapsed
res_collapsed <- NULL; for (i in 1:nrow(anu)) {
  tmp <- IBD_anubis[IBD_anubis$ID == anu$name[i] & IBD_anubis$source %in% anu_source$name[c(1:10,18)],]; tmp2 <- NULL
  for (chr in 1:20) {
    t3 <- tmp[tmp$chrom == chr,]; t3 <- t3[order(t3$start, t3$end),]
    cov <- coverage(IRanges(t3$start, t3$end)); as.data.table(cbind(cov@values, cov@lengths)) -> mc
    mc$end <- cumsum(mc$V2); mc$start <- c(1,mc$end[-length(mc$end)]+1); colnames(mc) <- c("number", "length", "end", "start"); mc$chrom <- chr
    tmp2 <- rbind(tmp2[,c(5,4,3,1,2)], mc); rm(mc,cov,t3)
  }; rm(chr)
  tmp2$name <- anu$name[i]; res_collapsed <- rbind(res_collapsed, tmp2); rm(tmp2)
}; rm(i, tmp, IBD_anubis); res_collapsed -> anu_collapsed; rm(res_collapsed)

# get just the segments which pass a threshold of 30, 50, or 70% of possible source individuals
# for yellow
ibd_mix_tracts <- NULL; min_n <- 8; for (i in 1:nrow(yel)) {
  tmp <- yel_collapsed[yel_collapsed$name == yel$name[i] & yel_collapsed$number > min_n,]
  tmp2 <- NULL
  for (chr in 1:20) {
    t3 <- tmp[tmp$chrom == chr,]
    t3 <- t3[order(t3$start, t3$end),]
    cov <- coverage(IRanges(t3$start, t3$end))
    as.data.table(cbind(cov@values, cov@lengths)) -> mc
    mc$end <- cumsum(mc$V2)
    mc$start <- c(1,mc$end[-length(mc$end)]+1)
    #mc <- subset(mc, mc$V1 > 0)
    colnames(mc) <- c("number", "length", "end", "start"); mc$chrom <- chr
    tmp2 <- rbind(tmp2, mc[,c(5,4,3,1,2)]); rm(mc,cov,t3)
  }; rm(chr); tmp2$name <- yel$name[i]
  ibd_mix_tracts <- rbind(ibd_mix_tracts, tmp2); rm(tmp2)
}; rm(i, tmp); ibd_mix_tracts -> ibd_mix_tracts_30yel; rm(ibd_mix_tracts)
ibd_mix_tracts <- NULL; min_n <- 14; for (i in 1:nrow(yel)) {
  tmp <- yel_collapsed[yel_collapsed$name == yel$name[i] & yel_collapsed$number > min_n,]
  tmp2 <- NULL
  for (chr in 1:20) {
    t3 <- tmp[tmp$chrom == chr,]
    t3 <- t3[order(t3$start, t3$end),]
    cov <- coverage(IRanges(t3$start, t3$end))
    as.data.table(cbind(cov@values, cov@lengths)) -> mc
    mc$end <- cumsum(mc$V2)
    mc$start <- c(1,mc$end[-length(mc$end)]+1)
    #mc <- subset(mc, mc$V1 > 0)
    colnames(mc) <- c("number", "length", "end", "start"); mc$chrom <- chr
    tmp2 <- rbind(tmp2, mc[,c(5,4,3,1,2)]); rm(mc,cov,t3)
  }; rm(chr); tmp2$name <- yel$name[i]
  ibd_mix_tracts <- rbind(ibd_mix_tracts, tmp2); rm(tmp2)
}; rm(i, tmp); ibd_mix_tracts -> ibd_mix_tracts_50yel; rm(ibd_mix_tracts)
ibd_mix_tracts <- NULL; min_n <- 19; for (i in 1:nrow(yel)) {
  tmp <- yel_collapsed[yel_collapsed$name == yel$name[i] & yel_collapsed$number > min_n,]
  tmp2 <- NULL
  for (chr in 1:20) {
    t3 <- tmp[tmp$chrom == chr,]
    t3 <- t3[order(t3$start, t3$end),]
    cov <- coverage(IRanges(t3$start, t3$end))
    as.data.table(cbind(cov@values, cov@lengths)) -> mc
    mc$end <- cumsum(mc$V2)
    mc$start <- c(1,mc$end[-length(mc$end)]+1)
    #mc <- subset(mc, mc$V1 > 0)
    colnames(mc) <- c("number", "length", "end", "start"); mc$chrom <- chr
    tmp2 <- rbind(tmp2, mc[,c(5,4,3,1,2)]); rm(mc,cov,t3)
  }; rm(chr); tmp2$name <- yel$name[i]
  ibd_mix_tracts <- rbind(ibd_mix_tracts, tmp2); rm(tmp2)
}; rm(i, tmp); ibd_mix_tracts -> ibd_mix_tracts_70yel; rm(ibd_mix_tracts)

# for anubis
ibd_mix_tracts <- NULL; min_n <- 3; for (i in 1:nrow(anu)) {
  tmp <- anu_collapsed[anu_collapsed$name == anu$name[i] & anu_collapsed$number > min_n,]
  tmp2 <- NULL
  for (chr in 1:20) {
    t3 <- tmp[tmp$chrom == chr,]
    t3 <- t3[order(t3$start, t3$end),]
    cov <- coverage(IRanges(t3$start, t3$end))
    as.data.table(cbind(cov@values, cov@lengths)) -> mc
    mc$end <- cumsum(mc$V2)
    mc$start <- c(1,mc$end[-length(mc$end)]+1)
    #mc <- subset(mc, mc$V1 > 0)
    colnames(mc) <- c("number", "length", "end", "start"); mc$chrom <- chr
    tmp2 <- rbind(tmp2, mc[,c(5,4,3,1,2)]); rm(mc,cov,t3)
  }; rm(chr); tmp2$name <- anu$name[i]
  ibd_mix_tracts <- rbind(ibd_mix_tracts, tmp2); rm(tmp2)
}; rm(i, tmp); ibd_mix_tracts -> ibd_mix_tracts_30anu; rm(ibd_mix_tracts)
ibd_mix_tracts <- NULL; min_n <- 5; for (i in 1:nrow(anu)) {
  tmp <- anu_collapsed[anu_collapsed$name == anu$name[i] & anu_collapsed$number > min_n,]
  tmp2 <- NULL
  for (chr in 1:20) {
    t3 <- tmp[tmp$chrom == chr,]
    t3 <- t3[order(t3$start, t3$end),]
    cov <- coverage(IRanges(t3$start, t3$end))
    as.data.table(cbind(cov@values, cov@lengths)) -> mc
    mc$end <- cumsum(mc$V2)
    mc$start <- c(1,mc$end[-length(mc$end)]+1)
    #mc <- subset(mc, mc$V1 > 0)
    colnames(mc) <- c("number", "length", "end", "start"); mc$chrom <- chr
    tmp2 <- rbind(tmp2, mc[,c(5,4,3,1,2)]); rm(mc,cov,t3)
  }; rm(chr)
  tmp2$name <- anu$name[i]
  ibd_mix_tracts <- rbind(ibd_mix_tracts, tmp2); rm(tmp2)
}; rm(i, tmp); ibd_mix_tracts -> ibd_mix_tracts_50anu; rm(ibd_mix_tracts)
ibd_mix_tracts <- NULL; min_n <- 7; for (i in 1:nrow(anu)) {
  tmp <- anu_collapsed[anu_collapsed$name == anu$name[i] & anu_collapsed$number > min_n,]
  tmp2 <- NULL
  for (chr in 1:20) {
    t3 <- tmp[tmp$chrom == chr,]
    t3 <- t3[order(t3$start, t3$end),]
    cov <- coverage(IRanges(t3$start, t3$end))
    as.data.table(cbind(cov@values, cov@lengths)) -> mc
    mc$end <- cumsum(mc$V2)
    mc$start <- c(1,mc$end[-length(mc$end)]+1)
    #mc <- subset(mc, mc$V1 > 0)
    colnames(mc) <- c("number", "length", "end", "start"); mc$chrom <- chr
    tmp2 <- rbind(tmp2, mc[,c(5,4,3,1,2)]); rm(mc,cov,t3)
  }; rm(chr)
  tmp2$name <- anu$name[i]
  ibd_mix_tracts <- rbind(ibd_mix_tracts, tmp2); rm(tmp2)
}; rm(i, tmp); ibd_mix_tracts -> ibd_mix_tracts_70anu; rm(ibd_mix_tracts)
rm(min_n)

save.image("./RESULTS/IBDmix_tracts.RData")

write.table(ibd_mix_tracts_30anu[ibd_mix_tracts_30anu$number > 0,], "./RESULTS/ibdmix_tracts.anubis.30.txt", row.names=F, col.names=T, sep="\t", quote=F)
write.table(ibd_mix_tracts_50anu[ibd_mix_tracts_50anu$number > 0,], "./RESULTS/ibdmix_tracts.anubis.50.txt", row.names=F, col.names=T, sep="\t", quote=F)
write.table(ibd_mix_tracts_70anu[ibd_mix_tracts_70anu$number > 0,], "./RESULTS/ibdmix_tracts.anubis.70.txt", row.names=F, col.names=T, sep="\t", quote=F)
write.table(ibd_mix_tracts_30yel[ibd_mix_tracts_30yel$number > 0,], "./RESULTS/ibdmix_tracts.yellow.30.txt", row.names=F, col.names=T, sep="\t", quote=F)
write.table(ibd_mix_tracts_50yel[ibd_mix_tracts_50yel$number > 0,], "./RESULTS/ibdmix_tracts.yellow.50.txt", row.names=F, col.names=T, sep="\t", quote=F)
write.table(ibd_mix_tracts_70yel[ibd_mix_tracts_70yel$number > 0,], "./RESULTS/ibdmix_tracts.yellow.70.txt", row.names=F, col.names=T, sep="\t", quote=F)
