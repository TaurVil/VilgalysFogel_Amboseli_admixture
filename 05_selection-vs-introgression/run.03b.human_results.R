## Get results for the human data which recapitulate those previously published
library(ppcor); library(DescTools)

load("./windows.human.250kb.RData")

## QC to see if thinning substantially alters results
## let's keep every other window
## all <- all[seq(2, nrow(all),2),]

## filter for windows which have Neanderthal genotypes for at least 60% of the sequence
## number of sites in the window was recorded as "testable"
all <- subset(all, all$testable > 150000)


## recombination results
t <- cor.test(all$mean_ancestry,all$recombination,method=c('spearman')); t; t$p.value;rm(t)
SpearmanRho(all$mean_ancestry,all$recombination, conf.level = 0.95)

## B value results
t <- cor.test(all$mean_ancestry,all$B,method=c('spearman')); t; t$p.value; rm(t)
SpearmanRho(all$mean_ancestry,all$B, conf.level = 0.95)

## Fixed differences
t <- cor.test(all$mean_ancestry,all$fixed_N3,method=c('spearman')); t; t$p.value; rm(t)
SpearmanRho(all$mean_ancestry,all$fixed_N3, conf.level = 0.95)

## Verify that results are robust to the number of SNPs which have Neanderthal genotypes
# library(ppcor)
# pcor(x=cbind(all$mean_ancestry, all$rcr, all$testable), method='spearman')$estimate
# pcor(x=cbind(all$mean_ancestry, all$B, all$testable), method='spearman')$estimate
# pcor(x=cbind(all$mean_ancestry, all$fixed_N3, all$testable), method='spearman')$estimate


## Recreate figure of fixed differences from Vernot and Akey

# a2 <- subset(all, all$mean_ancestry > 0)
# t <- cor.test(a2$mean_ancestry,a2$fixed,method=c('spearman')); t; t$p.value; rm(t)
# res2 <- as.data.frame(matrix(ncol=5, nrow=400))
# xxx <- 30
# for (i in seq(1,xxx+1)) {
#   if (i %in% 1:xxx) {z <- i-1; subset(a2, a2$fixed_N3 == z) -> tmp}
#   else {z <- i -1; subset(a2, a2$fixed_N3 >= z) -> tmp}
#   print(nrow(tmp))
#   
#   res2[i,1] <- z
#   res2[i,2] <- nrow(tmp)
#   res2[i,3] <- mean(tmp$mean_ancestry, na.rm=T)
#   
#   rm(tmp)
# }
# plot(res2$V3 ~ res2$V1, xlab="Number of fixed differences", ylab="mean Neanderthal ancestry", pch=16, cex=2, col='purple4'); cor.test(res2$V3, res2$V1, method='spearman')

