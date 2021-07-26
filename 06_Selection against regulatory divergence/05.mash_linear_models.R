## integrate local ancestry effects using adaptive shrinkage to have a single effect size estimate per gene
## we'll also look at the expression differences in each dataset, collapsing our previous permutations into a table of results with permuation-based FDRs

library(data.table); library(cobs); library(mashr); library(ggplot2); library(reshape2)

#######################################
## load in results
#######################################
dataset="eLife"; load(paste("./expression_results.", dataset, ".RData", sep=""))
dataset="TruCulture"; load(paste("./expression_results.", dataset, ".RData", sep=""))
rm(dataset)

#######################################
## filter for minimum standard error, get FDR
## sites with extremely low standard errors are produced by failures in model convergence and lead to erroneously significant p-values which don't appear to match the data
#######################################
min_se <- 10^-3.5
min_p <- 1e-50

eLife <- subset(eLife, eLife$se_local >= min_se & eLife$p_local >= min_p)
perm_eLife <- subset(perm_eLife, perm_eLife$se_local >= min_se & perm_eLife$p_local >= min_p)
eLife <- perm.fdr(eLife, perm_eLife$p_local, "p_local", "10perm")
genes_eLife <- genes_eLife[genes_eLife$gene_ID %in% eLife$gene,]

TC <- subset(TC, TC$se_local >= min_se & TC$p_local >= min_p)
perm_TC <- subset(perm_TC, perm_TC$se_local >= min_se & perm_TC$p_local >= min_p)
TC <- perm.fdr(TC, perm_TC$p_local, "p_local", "10perm")
genes_TC <- genes_TC[genes_TC$gene_ID %in% TC$gene,]

#######################################
## show permuted and observed distributions
## get number of genes which pass a 10% false discovery rate (0 for both datasets)
#######################################

par(mfrow=c(1,2)); hist(perm_eLife$p_local, col='lightblue', breaks=100, main="permuted eLife"); hist(perm_TC$p_local, col=rgb(1,0,0,.4), breaks=100, main="permuted TC");
par(mfrow=c(1,2))
t1 <- -log10(perm_eLife$p_local); t2 <- -log10(eLife$p_local)
qqplot(t1, t2, xlab="-log10(permutated)", ylab="-log10(observed)", xlim=c(0,10), ylim=c(0,10), main="eLife", frame=F); abline(a=0,b=1,col='red')
t1 <- -log10(perm_TC$p_local); t2 <- -log10(TC$p_local)
qqplot(t1, t2, xlab="-log10(permutated)", ylab="-log10(observed)", xlim=c(0,10), ylim=c(0,10), main="TruCulture", frame=F); abline(a=0,b=1,col='red')

sum(TC$fdr_10perm < 0.1) ## 1,750
sum(eLife$fdr_10perm < 0.1) ## 878

#######################################
## filter for genes analyzed in each dataset
## 10,222 genes analyzed in both datasets. 
#######################################
genes <- genes_TC[genes_TC$gene_ID %in% genes_eLife$gene_ID,]
genes_TC <- genes_TC[!genes_TC$gene_ID %in% genes$gene_ID,]
genes_eLife <- genes_eLife[!genes_eLife$gene_ID %in% genes$gene_ID,]

only_eLife <- subset(eLife, !eLife$gene %in% genes$gene_ID)
only_TC <- subset(TC, !TC$gene %in% genes$gene_ID)

eLife <- subset(eLife, eLife$gene %in% genes$gene_ID)
TC <- subset(TC, TC$gene %in% genes$gene_ID)

eLife <- eLife[order(eLife$gene),]
TC <- TC[order(TC$gene),]
genes <- genes[order(genes$gene_ID),]

nrow(eLife); nrow(TC) ## tested 10,212 genes
sum(eLife$fdr_10perm < 0.1 | TC$fdr_10perm < 0.1) ## 2,047 genes

#######################################
## format combined dataset which has Betas and the SE_betas
#######################################
for_mash <- NULL; for_mash$Bhat <- for_mash$Shat <- matrix(ncol=2,nrow=nrow(genes))
row.names(for_mash$Bhat) <- row.names(for_mash$Shat) <- genes$Gene
colnames(for_mash$Bhat) <- colnames(for_mash$Shat) <- c("eLife", "TC")
for_mash$Bhat[,1] <- eLife$beta_local; for_mash$Bhat[,2] <- TC$beta_local
for_mash$Shat[,1] <- sqrt(eLife$se_local); for_mash$Shat[,2] <- sqrt(TC$se_local)

data=mash_set_data(for_mash$Bhat, for_mash$Shat)

#######################################
## plot correlation in effect sizes 
## they are highly correlated, consistent with a true genetic effect
#######################################
par(mfrow=c(1,1))
plot(for_mash$Bhat[,1] ~ for_mash$Bhat[,2], pch=16, col=rgb(0,0,1,.15), 
     frame=F, ylab="beta(eLife)", xlab="beta(TruCulture)")
abline(lm(for_mash$Bhat[,1] ~ for_mash$Bhat[,2]), col='red')
model <- summary(lm(for_mash$Bhat[,1] ~ for_mash$Bhat[,2])); model
legend("topleft", legend = paste("r2 = ", signif(model$r.squared, digits=2), "\nslope = ", signif(model$coefficients[2,1], digits = 2),"\neLife has effect sizes weaker than TC"), bty='n')
rm(model)

#######################################
## create possible covariances matrices, and run mash to adjust Bs
#######################################
U.c = cov_canonical(data) ## canonical covariance
## data driven covariance
m.1by1 = mash_1by1(data)
strong = get_significant_results(m.1by1,0.05)
U.pca = cov_pca(data,2,subset=strong)
U.ed = cov_ed(data, U.pca, subset=strong)

# run mash
m  = mash(data, c(U.c,U.ed,U.pca))


#######################################
## print which covariance matrices dominate the expression data, and a histogram of local false sign rates after mash
#######################################
tmp <- as.data.frame(get_estimated_pi(m)); colnames(tmp) <- 'v1'
row.names(tmp) -> tmp$name
ggplot(tmp, aes(x=name, y=v1)) + 
  geom_bar(stat="identity") + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle=30,hjust=1,vjust=1.0))

par(mfrow=c(1,2))
hist(get_lfsr(m)[,1], xlab="local false sign rate", main="eLife")
hist(get_lfsr(m)[,2], xlab="local false sign rate", main="TC")
dev.off()

#######################################
## get effect sizes from mash
#######################################
cbind(get_pm(m), get_psd(m), get_lfsr(m)) -> mash_res
colnames(mash_res)[3:4] <- c("SE_eLife", "SE_TC")
colnames(mash_res)[5:6] <- c("lfsr_eLife", "lfsr_TC")
as.data.frame(mash_res) -> mash_res
# mash_res$abs_eLife <- abs(mash_res$eLife)
# mash_res$abs_TC <- abs(mash_res$TC)
mash_res$gene <- genes$gene_ID
mash_res$abs_mean <- abs((mash_res$eLife + mash_res$TC)/2)
# mash_res$max_abs <- apply(mash_res[,5:6],1,max)

TC$fdr_10perm -> mash_res$q_TC
eLife$fdr_10perm -> mash_res$q_eLife
eLife$beta_local -> mash_res$eLife_rawB
TC$beta_local -> mash_res$TC_rawB
eLife$se_local -> mash_res$eLife_rawSE
TC$se_local -> mash_res$TC_rawSE
eLife$p_local -> mash_res$eLife_rawP
TC$p_local -> mash_res$TC_rawP

#######################################
## save table of results
#######################################
write.table(mash_res, "./mash_results.txt", row.names=F, col.names=T, quote=F, sep="\t")

## Local false sign rates estimated from mash are largely concordant with our FDR-based estimates
sum(mash_res$lfsr_eLife < 0.1 | mash_res$lfsr_TC < 0.1); nrow(mash_res)
sum(mash_res$lfsr_eLife < 0.1 | mash_res$lfsr_TC < 0.1)/nrow(mash_res)

## Correlation in effect size after mash
summary(lm(mash_res$eLife_rawB ~ mash_res$TC_rawB))
summary(lm(mash_res$eLife_rawB ~ mash_res$TC_rawB))$coefficients[2,4]

## Correlation in effect size before mash
summary(lm(mash_res$eLife ~ mash_res$TC))
summary(lm(mash_res$eLife ~ mash_res$TC))$coefficients[2,4]

