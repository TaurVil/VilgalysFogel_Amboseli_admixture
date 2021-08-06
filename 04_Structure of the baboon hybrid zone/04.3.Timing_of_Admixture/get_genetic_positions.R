#!/usr/bin/env Rscript
#SBATCH --get-user-env

## read in snp file without position in cM
library(data.table); fread("ref.snp") -> snp
s2 <- as.data.frame(cbind(paste("snp", snp$V4,".", snp$V2, snp$V5, snp$V6,sep=""), snp$V2))
colnames(s2)[1:2] <- c("snpid","chr")
s2$GP <- NA
s2$PP <- snp$V4

## load in recombination rates. This file is produced in 05_selection_against_introgression. 
load("./Windows_25kb_recombination.RData")
## removes windows with unreasonable recombination
subset(rcr, rcr$n24_anubis < 100*median(rcr$n24_anubis)) -> res2
s3 <- NULL 

## for each chromosome, get SNPs, recombination rate, and total length. We'll then fit a smoothed function to avoid big jumps in the genetic distance which are implausible and produce a vector of genetic positions to accompany the physical positions. 
for (i in unique(snp$V2)) {
	tmp=subset(s2, s2$chr == i) 
	r1=subset(res2, res2$chr == paste("chr",i,sep=""))
	r1$sum <- cumsum(r1$n24_anubis)
	
	d <- c(172,115,100,164,113,111,116,93,88,88,91,126,77,76,69,89,77,84,63,63) [i]
	
	r1$cM <- r1$sum*d/max(r1$sum)
	lo <- loess(r1$cM ~ r1$end, span = 0.01)
	tmp$GP <- predict(lo, tmp$PP)
	tmp$GP[tmp$PP < 75000] <- predict(lo,75000)*tmp$PP[tmp$PP < 75000]/ 75000
	tmp <- tmp[!is.na(tmp$GP),]
	rbind(s3,tmp) -> s3
	print(i)
}
#s2$GP <- s2$PP/1e6
write.table(s3, "n46.v2.snp", row.names=F, col.names=F, sep="\t", quote=F)
