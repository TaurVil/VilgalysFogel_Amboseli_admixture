#!/usr/bin/env Rscript
#SBATCH --get-user-env

library(data.table); fread("family_packed.snp") -> snp
fread("ambo.v2.snp") -> tst
tst$V3/100 -> tst$V3 # switch from cM to M 
subset(tst, tst$V1 %in% snp$V1) -> t2 # only include SNPs called in Amboseli
t2$V3 -> snp$V3; rm(t2,tst)

snp$V1 <- paste("snp",row.names(snp),sep="")
s2 <- NULL 
for (i in unique(snp$V2)) {
	s=subset(snp, snp$V2 == i)
	s<-s[order(s$V4),]
	s$p <- c(-0.00025,s$V3[-nrow(s)])
	s$V3-s$p -> s$d 
	s$d[s$d < 1e-16] <- 1e-16
	s$V3 <- cumsum(s$d)
	#sum(order(s$V3)==order(s$V4))/nrow(s)
	rbind(s2,s) -> s2; print(i) 
}; write.table(s2[,1:6], "family_packed.snp", row.names=F, col.names=F, sep="\t", quote=F)


