#!/usr/bin/env Rscript
#SBATCH --get-user-env

library(data.table); library(plyr)
fread("./true_tracts/NAME.tracts.txt") -> sample
read.delim("vcf_example_header") -> VCF

subset(sample, sample$chrom == "SCAF") -> sample

nom=paste("SCAF",".anubis.frq",sep="")
anu <- fread(nom); rm(nom); anu$V1 <- paste("chr",anu$V1,sep="")
colnames(anu) <- c('chr', 'pos','n_alleles','n_chr','ref','ref_freq','alt','alt_freq')
nom=paste("SCAF",".yellow.frq",sep="")
yel <- fread(nom); rm(nom); yel$V1 <- paste("chr",yel$V1,sep="")
colnames(yel) <- c('chr', 'pos','n_alleles','n_chr','ref','ref_freq','alt','alt_freq')

subset(anu, !(anu$alt == "*")) -> anu
subset(yel, !(yel$alt == "*")) -> yel

vcf <- as.data.frame(matrix(ncol=10, nrow=nrow(anu))); colnames(vcf) <- colnames(VCF)
k=1; unique(c(anu$pos, yel$pos)) -> sites

set.seed(100)

colnames(vcf)[1] <- "CHROM"
vcf$CHROM <- "SCAF"
vcf$POS <- sites
vcf$REF <- anu$ref
vcf$ALT <- anu$alt
vcf$INFO <- vcf$ID <- "."
vcf$QUAL <- 100
vcf$FILTER <- "QD"
vcf$FORMAT <- "GT"


get_gt <- function(site) {
        ancestry <- subset(sample, sample$start <= site & sample$end > site)$ancestry
        anu_freq <- anu[anu$pos == site,]$alt_freq
        yel_freq <- yel[yel$pos == site,]$alt_freq

        if (!is.na(anu_freq) & !is.na(yel_freq) & length(ancestry) > 0) {
                # make two calls from yellow
                if (ancestry == 0) {
                        jk <- sum(runif(2,0,1) < yel_freq)
                }
                # make two calls from anubis
                if (ancestry == 2) {
                        jk <- sum(runif(2,0,1) < anu_freq)
                }
                # make one call from each
                if (ancestry == 1) {
                        jk <- sum(runif(1,0,1) < yel_freq, runif(1,0,1) < anu_freq)
                }

                if (jk == 0) {genotype <- "0/0"}; if (jk == 1) {genotype <- "0/1"}; if (jk == 2) {genotype <- "1/1"}

        } else {genotype <- NA}

        return(genotype)
}

library(parallel)
vcf[,10] <- do.call("c",mclapply(sites,get_gt))
colnames(vcf)[10] <- "NAME"


subset(vcf, !is.na(vcf[,10])) -> vcf

write.table(vcf, "./simulated_vcfs/NAME.SCAF.vcf", row.names=F, col.names=T, quote=F, sep="\t")
