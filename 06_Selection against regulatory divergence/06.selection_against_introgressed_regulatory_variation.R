## hypothesis: selection is strongest against genes with increased regulatory divergence
#######################################

library(data.table); library(ggplot2); library(parallel); library(DescTools)

#######################################
## load in expression results (after mash) and recombination rate/local ancestry per gene
#######################################

# expression results
expr_results <- read.delim("./mash_results.txt")

# masked vs unmasked refpanel 
load("./Genes_ancestry.RData")

#######################################
## Recapitulate recombination-rate ancestry results for all genes
#######################################
t <- cor.test(features_genes$mean_ancestry,features_genes$recombination_rate,method=c('spearman')); t; t$p.value; rm(t)
t <- cor.test(features_genes$historical_ancestry,features_genes$recombination_rate,method=c('spearman')); t; t$p.value; rm(t)
t <- cor.test(features_genes$recent_ancestry,features_genes$recombination_rate,method=c('spearman')); t; t$p.value; rm(t)
######

#######################################
## Get just the genes which we analyzed expression data for (10,222 genes)
#######################################
to_analyze <- subset(features_genes, features_genes$gene_id %in% expr_results$gene); rm(features_genes)
expr_results <- subset(expr_results, expr_results$gene %in% to_analyze$gene_id)
to_analyze <- to_analyze[order(to_analyze$gene_id),]
expr_results <- expr_results[order(expr_results$gene),]

sum(to_analyze$gene_id == expr_results$gene) ## 10,222

#######################################
## Recapitulate recombination-rate ancestry results for analyzed genes
#######################################
t <- cor.test(to_analyze$mean_ancestry,to_analyze$recombination_rate,method=c('spearman')); t; t$p.value; rm(t)
t <- cor.test(to_analyze$historical_ancestry,to_analyze$recombination_rate,method=c('spearman')); t; t$p.value; rm(t)
t <- cor.test(to_analyze$recent_ancestry,to_analyze$recombination_rate,method=c('spearman')); t; t$p.value; rm(t)
######

#######################################
## Genes with large ancestry effects have lower levels of anubis ancestry
#######################################

## get ancestry values (already loaded in as "ancestry_genes") for just the genes where expression was compared
ancestry_genes <- ancestry_genes[ancestry_genes$gene_id %in% expr_results$gene,]
ancestry_genes <- ancestry_genes[order(ancestry_genes$gene_id),]
sum(ancestry_genes$gene_id == expr_results$gene) ## 10,222

## for the most extreme X% of genes (quantile), what is the difference in anubis ancestry across individuals? (paired t-test)
res_indiv <- as.data.frame(matrix(ncol=9, nrow=7))
colnames(res_indiv) <- c("quantile", "p", "diff", "CI_min", "CI_max", 
                         "p_historic", "diff_historic" ,"his_min", "his_max")
z<- 1; for (qqq in seq(0.1,0.4,0.05)) {
  res_indiv$quantile[z] <- qqq
  
  most <- subset(ancestry_genes, expr_results$abs_mean >= quantile(expr_results$abs_mean, 1-qqq))
  least <- subset(ancestry_genes, expr_results$abs_mean <= quantile(expr_results$abs_mean, qqq))
  
  model <- t.test(colMeans(most[,-c(1:5)]), colMeans(least[,-c(1:5)]), paired = T)
  model$p.value -> res_indiv$p[z]; model$estimate -> res_indiv$diff[z]; model$conf.int[1:2] -> res_indiv[z,4:5]
  
  hist_ids <- ids$table_s1_id[ids$recent_hybridsc == "historical"]; model <- t.test(colMeans(most[,colnames(most) %in% hist_ids]), colMeans(least[,colnames(least) %in% hist_ids]), paired = T)
  model$p.value -> res_indiv$p_historic[z]
  model$estimate -> res_indiv$diff_historic[z]
  model$conf.int[1:2] -> res_indiv[z,8:9]
  
  z <- z+1
}; rm(z, qqq, model, least, most, hist_ids, genes, ancestry_genes)

res_indiv
write.table(res_indiv, "./reduced_ancestry_in_DE_genes.txt", row.names=F, col.names=T, sep="\t", quote=F)

#######################################
## Genes with large ancestry effects have stronger signatures of selection (a stronger correlation with recombination rate)
## Calculate rho between ancestry and recombination for top/bottom quantiles
#######################################
res <- as.data.frame(matrix(ncol=9, nrow=7))
colnames(res) <- c("quantile", "rho_DE", "low_DE", "hi_DE", "rho_nDE", "low_nDE", "hi_nDE", "pval_DE", "pval_nDE")
z<- 1; for (qqq in seq(0.1,0.4,0.05)) {
  res$quantile[z] <- qqq
  DE <- subset(expr_results, expr_results$abs_mean >= quantile(expr_results$abs_mean,1-qqq))
  subset(to_analyze, to_analyze$gene_id %in% DE$gene) -> DE
  res[z,2:4] <- SpearmanRho(DE$historical_ancestry, DE$recombination_rate, conf.level = 0.95)
  res[z,8] <- cor.test(DE$historical_ancestry, DE$recombination_rate, method="spearman")$p.value
  
  nDE <- subset(expr_results, expr_results$abs_mean <= quantile(expr_results$abs_mean,qqq))
  subset(to_analyze, to_analyze$gene_id %in% nDE$gene) -> nDE
  res[z,5:7] <- SpearmanRho(nDE$historical_ancestry, nDE$recombination_rate, conf.level = 0.95)
  res[z,9] <- cor.test(nDE$historical_ancestry, nDE$recombination_rate, method="spearman")$p.value
  
  z <- z+1
}; rm(z, qqq, nDE, DE)

# Bootstrap top/bottom DE genes and compare rho
# qqq sets quantile, if we're not going through a loop
res$bootstrap_pval <- res$sd_DE <- res$sd_nDE <- NA
for (qqq in res$quantile) {
  print(paste("Quantile: ", qqq*100,"%", sep=""))
  DE <- subset(expr_results, expr_results$abs_mean >= quantile(expr_results$abs_mean,1-qqq))
  subset(to_analyze, to_analyze$gene_id %in% DE$gene) -> DE
  
  nDE <- subset(expr_results, expr_results$abs_mean <= quantile(expr_results$abs_mean,qqq))
  subset(to_analyze, to_analyze$gene_id %in% nDE$gene) -> nDE
  
  set.seed(42)
  boot_DE <- NULL; for (i in 1:10000) {
    tmp_DE <- DE[sample(1:nrow(DE), nrow(DE), replace=T),]
    tmp_rho <- summary(lm(rank(tmp_DE$historical_ancestry)~rank(tmp_DE$recombination_rate)))
    c(boot_DE, tmp_rho$coefficients[2,1]) -> boot_DE
  }; rm(tmp_DE, tmp_rho, i)
  
  boot_nDE <- NULL; for (i in 1:10000) {
    tmp_nDE <- nDE[sample(1:nrow(nDE), nrow(nDE), replace=T),]
    tmp_rho <- summary(lm(rank(tmp_nDE$historical_ancestry)~rank(tmp_nDE$recombination_rate)))
    c(boot_nDE, tmp_rho$coefficients[2,1]) -> boot_nDE
  }; rm(tmp_nDE, tmp_rho, i)
  
  sd(boot_DE) -> res$sd_DE[res$quantile == qqq]
  sd(boot_nDE) -> res$sd_nDE[res$quantile == qqq]
  
  diff <- boot_DE - boot_nDE
  print(qqq, quantile(diff, c(0.05,0.95))); print(ecdf(diff)(0))
  ecdf(diff)(0) -> res$bootstrap_pval[res$quantile == qqq]
}; rm(DE, nDE, qqq, diff, boot_DE, boot_nDE, expr_results)
res
write.table(res, "./bootstrap_rho_results.txt", row.names=F, col.names=T, sep="\t", quote=F)


