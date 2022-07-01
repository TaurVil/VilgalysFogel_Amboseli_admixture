## Ancestry and recombination rates 
library(data.table); library(ggplot2); library(parallel); library(DescTools); library(ppcor)

#########   Chunk for metadata       ###########

#########   Chunk to get data       ###########
load("../VilgalysFogel_main_data_file.250kb_windows.RData")

#########   Recombination rate result       ###########
t <- cor.test(to_analyze$mean_ancestry,to_analyze$rcr,method=c('spearman')); t; t$p.value; rm(t)
t <- cor.test(to_analyze$recent_ancestry,to_analyze$rcr,method=c('spearman')); t; t$p.value; rm(t)
t <- cor.test(to_analyze$historical_ancestry,to_analyze$rcr,method=c('spearman')); t; t$p.value; rm(t)

#########   B values       ###########
t <- cor.test(to_analyze$mean_ancestry,to_analyze$B,method=c('spearman')); t; t$p.value; rm(t)
t <- cor.test(to_analyze$recent_ancestry,to_analyze$B,method=c('spearman')); t; t$p.value; rm(t)
t <- cor.test(to_analyze$historical_ancestry,to_analyze$B,method=c('spearman')); t; t$p.value; rm(t)

#########   number of variable sites      ###########
cor.test(to_analyze$mean_ancestry, to_analyze$snps, method="spearman")

#########   number of fixed differences      ###########
t <- cor.test(to_analyze$mean_ancestry, to_analyze$fixed, method="spearman"); t; t$p.value
t <- cor.test(to_analyze$recent_ancestry, to_analyze$fixed, method="spearman"); t; t$p.value
t <- cor.test(to_analyze$historical_ancestry, to_analyze$fixed, method="spearman"); t; t$p.value

#########  Fst, varying threshold (Supplementary Figure using 250 kb) ###########
results <- as.data.frame(matrix(ncol=6, nrow=11)) ; 
colnames(results) <- c("fixed_differences", "rho", "lower", "upper", "n_snps", "pvalue")
results[,1] <- seq(0.5,1,0.05)

{ # rho
  results[11,2:4] <- SpearmanRho(to_analyze$mean_ancestry, to_analyze$fixed, conf.level = 0.95)
  results[10,2:4] <- SpearmanRho(to_analyze$mean_ancestry, to_analyze$fst95, conf.level = 0.95)
  results[9,2:4] <- SpearmanRho(to_analyze$mean_ancestry, to_analyze$fst90, conf.level = 0.95)
  results[8,2:4] <- SpearmanRho(to_analyze$mean_ancestry, to_analyze$fst85, conf.level = 0.95)
  results[7,2:4] <- SpearmanRho(to_analyze$mean_ancestry, to_analyze$fst80, conf.level = 0.95)
  results[6,2:4] <- SpearmanRho(to_analyze$mean_ancestry, to_analyze$fst75, conf.level = 0.95)
  results[5,2:4] <- SpearmanRho(to_analyze$mean_ancestry, to_analyze$fst70, conf.level = 0.95)
  results[4,2:4] <- SpearmanRho(to_analyze$mean_ancestry, to_analyze$fst65, conf.level = 0.95)
  results[3,2:4] <- SpearmanRho(to_analyze$mean_ancestry, to_analyze$fst60, conf.level = 0.95)
  results[2,2:4] <- SpearmanRho(to_analyze$mean_ancestry, to_analyze$fst55, conf.level = 0.95)
  results[1,2:4] <- SpearmanRho(to_analyze$mean_ancestry, to_analyze$fst50, conf.level = 0.95)
}
{ # pvalue
  results[11,6] <- cor.test(to_analyze$mean_ancestry, to_analyze$fixed, method='spearman')$p.value
  results[10,6] <- cor.test(to_analyze$mean_ancestry, to_analyze$fst95, method='spearman')$p.value
  results[9,6] <- cor.test(to_analyze$mean_ancestry, to_analyze$fst90, method='spearman')$p.value
  results[8,6] <- cor.test(to_analyze$mean_ancestry, to_analyze$fst85, method='spearman')$p.value
  results[7,6] <- cor.test(to_analyze$mean_ancestry, to_analyze$fst80, method='spearman')$p.value
  results[6,6] <- cor.test(to_analyze$mean_ancestry, to_analyze$fst75, method='spearman')$p.value
  results[5,6] <- cor.test(to_analyze$mean_ancestry, to_analyze$fst70, method='spearman')$p.value
  results[4,6] <- cor.test(to_analyze$mean_ancestry, to_analyze$fst65, method='spearman')$p.value
  results[3,6] <- cor.test(to_analyze$mean_ancestry, to_analyze$fst60, method='spearman')$p.value
  results[2,6] <- cor.test(to_analyze$mean_ancestry, to_analyze$fst55, method='spearman')$p.value
  results[1,6] <- cor.test(to_analyze$mean_ancestry, to_analyze$fst50, method='spearman')$p.value
}

ggplot(data=results, aes(x=fixed_differences, y=rho)) + 
  geom_segment(aes(x=fixed_differences, xend=fixed_differences, y=lower, yend=upper), color = "darkgreen", size=1.4) +
  geom_point(col = 'darkgreen', size=8) + 
  theme_classic() + xlab("FST threshold") + ylab("Spearman's rho") + #ylim(c(-0.151,0)) + 
  theme(axis.text = element_text(color="black"), text=element_text(size=16))

rm(t, results)
