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
min_se <- 10^-2
min_p <- 1e-8

eLife <- subset(eLife, eLife$se_global >= min_se & eLife$p_global >= min_p)
perm_global_eLife <- subset(perm_global_eLife, perm_global_eLife$se_global >= min_se & perm_global_eLife$p_global >= min_p)
eLife <- perm.fdr(eLife, perm_global_eLife$p_global, "p_global", "10perm")
genes_eLife <- genes_eLife[genes_eLife$gene_ID %in% eLife$gene,]

TC <- subset(TC, TC$se_global >= min_se & TC$p_global >= min_p)
perm_global_TC <- subset(perm_global_TC, perm_global_TC$se_global >= min_se & perm_global_TC$p_global >= min_p)
TC <- perm.fdr(TC, perm_global_TC$p_global, "p_global", "10perm")
genes_TC <- genes_TC[genes_TC$gene_ID %in% TC$gene,]

#######################################
## show permuted and observed distributions
## get number of genes which pass a 10% false discovery rate (0 for both datasets)
#######################################

## null distributions are slightly enriched for apparent effects of global ancestry, perhaps a consequence of cryptic structure
par(mfrow=c(1,2)); hist(perm_global_eLife$p_global, col='lightblue', breaks=100, main="permuted eLife"); hist(perm_global_TC$p_global, col=rgb(1,0,0,.4), breaks=100, main="permuted TC");

## observed results for global ancestry show no enrichment of significant effects
par(mfrow=c(1,2))
t1 <- -log10(perm_global_eLife$p_global); t2 <- -log10(eLife$p_global)
qqplot(t1, t2, xlab="-log10(permutated)", ylab="-log10(observed)", main="eLife", frame=F); abline(a=0,b=1,col='red')
t1 <- -log10(perm_global_TC$p_global); t2 <- -log10(TC$p_global)
qqplot(t1, t2, xlab="-log10(permutated)", ylab="-log10(observed)", main="TruCulture", frame=F); abline(a=0,b=1,col='red')
rm(t1, t2, perm.fdr, perm_global_TC, perm_global_eLife, perm_TC, perm_eLife)

sum(TC$fdr_10perm < 0.1)
sum(eLife$fdr_10perm < 0.1)

#######################################
## filter for genes analyzed in each dataset
## 10,222 genes analyzed in both datasets. None have significant effects of global ancestry. 
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

sum(eLife$fdr_10perm < 0.1 | TC$fdr_10perm < 0.1)

#######################################
## format combined dataset which has Betas and the SE_betas
#######################################
for_mash <- NULL; for_mash$Bhat <- for_mash$Shat <- matrix(ncol=2,nrow=nrow(genes))
row.names(for_mash$Bhat) <- row.names(for_mash$Shat) <- genes$Gene
colnames(for_mash$Bhat) <- colnames(for_mash$Shat) <- c("eLife", "TC")
for_mash$Bhat[,1] <- eLife$beta_global; for_mash$Bhat[,2] <- TC$beta_global
for_mash$Shat[,1] <- sqrt(eLife$se_global); for_mash$Shat[,2] <- sqrt(TC$se_global)

data=mash_set_data(for_mash$Bhat, for_mash$Shat)

#######################################
## plot correlation in effect sizes 
## they're not strongly correlated, which is consistent with no real biological signal
#######################################
par(mfrow=c(1,1))
plot(for_mash$Bhat[,1] ~ for_mash$Bhat[,2], pch=16, col=rgb(0,0,1,.15), 
     frame=F, ylab="beta(eLife)", xlab="beta(TruCulture)")
abline(lm(for_mash$Bhat[,1] ~ for_mash$Bhat[,2]), col='red')



