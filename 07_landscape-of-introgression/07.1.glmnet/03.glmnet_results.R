## Read in output file from 02
load("./03.glmnet_model_results.R")

## predictive power of population-wide ancestry 
mean(sqrt(results$r2_stand)); sd(sqrt(results$r2_stand))
## predictive power of permutations
mean(sqrt(results_perm$r2_stand)); sd(sqrt(results_perm$r2_stand))

## predictive power of historical anubis ancestry
mean(sqrt(results_hist$r2_stand)); sd(sqrt(results_hist$r2_stand))
## predictive power of recent anubis ancestry
mean(sqrt(results_recent$r2_stand)); sd(sqrt(results_recent$r2_stand))

## convert commonly used predictors into categories for Fig 4C
b2 <- betas[, colSums(betas == 0)/nrow(betas) < 1]
as.data.frame(b2[,-1]) -> b2
head(b2)
b2$`background selection` <- rowSums(b2[,colnames(b2) %in% c("B_exons", "B")])
b2 <- b2[,-c(1,7)]
b2$`recombination rate` <- rowSums(b2[,colnames(b2) %in% c("rankrcr", "logrcr")])
b2 <- b2[,-c(6,7)]
colnames(b2)[colnames(b2) == 'snps'] <- "number of SNVs"
b2$`highly differentiated sites` <- rowSums(b2[,2:5])
b2 <- b2[,-c(2:5)]

write.table(b2, "./glmnet_results.txt", row.names=T, col.names=T, sep="\t", quote=F)
