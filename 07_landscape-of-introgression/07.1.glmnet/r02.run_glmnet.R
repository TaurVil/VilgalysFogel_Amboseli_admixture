###############################
## Plot covariance
###############################
load("./ancestry_and_features.RData")

library(ggplot2)
matcov <- cor(data, method = "spearman"); library(reshape2); 
# matcov[,] -> matcov # reorder for a cleaner plot? 
melt(matcov) -> melt_matcov
ggplot(data = melt_matcov, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + 
  theme(axis.text.x = element_text(angle=90)) + 
  labs(x="", y="")
matcov <- cor(data, method = "pearson"); library(reshape2); 
melt(matcov) -> melt_matcov
ggplot(data = melt_matcov, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + 
  theme(axis.text.x = element_text(angle=90)) + 
  labs(x="", y="")
rm(matcov, melt_matcov)

###############################
## Run glmnet, mean population-wide ancestry 
###############################
library(glmnet)
alpha = 0.5 # elastic net parameter. the model's predictive power was optimized when alpha = 0.5, analyses not included. 

# Define response variable and predictors
x <- as.matrix(data[,-c(2:4)]); y <- as.matrix(data$mean_ancestry)

results <- betas <- NULL
set.seed(42); for (rep in 1:200) {
  ## Get set of regions to model (75%) and to test (25%)
  sample(x = 1:nrow(data), size = nrow(data)*0.75) -> list_test 
  test_x <- as.matrix(x[list_test,])
  test_y <- as.matrix(y[list_test,]); log(test_y) -> test_y
  pred_x <- as.matrix(x[-list_test,])
  pred_y <- as.matrix(y[-list_test,]); log(pred_y) -> pred_y
  
  print(rep)

    ## create data frame to index results
    tmp_results <- as.data.frame(matrix(ncol=5))
    colnames(tmp_results) <- c("rep", "alpha", "lambda", "r2_stand")
    
    tmp_results$rep <- rep; tmp_results$alpha <- alpha
    
    ## run glmnet
    model<-cv.glmnet(test_x,test_y,nfolds=20,alpha=alpha,standardize=T)
    tmp_results$lambda <- model$lambda.min
    ## predict ancestry of the test_set
    predicted <- predict(model, newx=pred_x, s="lambda.min")
    tmp_results$r2_stand <- summary(lm(predicted ~ pred_y))$r.squared
    weights<-t(as.matrix(unlist(coef(model,lambda="lambda.min"))[,1], nrow=1, ncol=30))
    rm(model, predicted)
    
    ## merge results with previous iterations
    rbind(results, tmp_results) -> results; rm(tmp_results)
    rbind(betas, weights) -> betas; rm(weights)
   
  # clean the working environment before exiting
  rm(alpha, pred_x, pred_y, test_x, test_y)
  print("done")
}; rm(rep)

results_perm <- NULL
set.seed(42); for (rep in 1:200) {
  ## Get set of regions to model and to test
  sample(x = 1:nrow(data), size = nrow(data)*0.75) -> list_test
  matrix(ncol=1, sample(unlist(y))) -> new_y
  test_x <- as.matrix(x[list_test,])
  test_y <- as.matrix(new_y[list_test,]); log(test_y) -> test_y
  pred_x <- as.matrix(x[-list_test,])
  pred_y <- as.matrix(new_y[-list_test,]); log(pred_y) -> pred_y
  
  print(rep)
  
  alpha=0.5
  tmp_results <- as.data.frame(matrix(ncol=5))
  colnames(tmp_results) <- c("rep", "alpha", "lambda", "r2_stand")
  
  tmp_results$rep <- rep; tmp_results$alpha <- alpha
  
  model<-cv.glmnet(test_x,test_y,nfolds=20,alpha=alpha,standardize=T)
  tmp_results$lambda <- model$lambda.min
  predicted <- predict(model, newx=pred_x, s="lambda.min")
  tmp_results$r2_stand <- summary(lm(predicted ~ pred_y))$r.squared
  
  rm(model, predicted)
  rbind(results_perm, tmp_results) -> results_perm; rm(tmp_results)
  
  rm(alpha)
  rm(pred_x, pred_y, test_x, test_y, new_y)
  print("done")
}; rm(rep)
## remove permuted results with no features that are unable to converge and erroneously result in near perfect fits
subset(results_perm, results_perm$r2_stand < 0.4) -> results_perm # we need to get rid of sets with no features because the models are unable to converge

###############################
## Run glmnet, recent and historical ancestry 
###############################

y <- as.matrix(data$recent_ancestry)

results_recent <- betas_recent <- NULL
set.seed(42); for (rep in 1:200) {
  ## Get set of regions to model and to test
  sample(x = 1:nrow(data), size = nrow(data)*0.75) -> list_test 
  test_x <- as.matrix(x[list_test,])
  test_y <- as.matrix(y[list_test,]); log(test_y) -> test_y
  pred_x <- as.matrix(x[-list_test,])
  pred_y <- as.matrix(y[-list_test,]); log(pred_y) -> pred_y
  
  print(rep)
  
  alpha=0.5
  tmp_results <- as.data.frame(matrix(ncol=5))
  colnames(tmp_results) <- c("rep", "alpha", "lambda", "r2_stand")
  
  tmp_results$rep <- rep; tmp_results$alpha <- alpha
  
  model<-cv.glmnet(test_x,test_y,nfolds=20,alpha=alpha,standardize=T)
  tmp_results$lambda <- model$lambda.min
  predicted <- predict(model, newx=pred_x, s="lambda.min")
  tmp_results$r2_stand <- summary(lm(predicted ~ pred_y))$r.squared
  weights<-t(as.matrix(unlist(coef(model,lambda="lambda.min"))[,1], nrow=1, ncol=30))
  
  rm(model, predicted)
  rbind(results_recent, tmp_results) -> results_recent; rm(tmp_results)
  rbind(betas_recent, weights) -> betas_recent; rm(weights)
  
  rm(alpha)
  rm(pred_x, pred_y, test_x, test_y)
  print("done")
}; rm(rep)

y <- as.matrix(data$historical_ancestry)

results_hist <- betas_hist <- NULL
set.seed(42); for (rep in 1:200) {
  ## Get set of regions to model and to test
  sample(x = 1:nrow(data), size = nrow(data)*0.75) -> list_test 
  test_x <- as.matrix(x[list_test,])
  test_y <- as.matrix(y[list_test,]); log(test_y) -> test_y
  pred_x <- as.matrix(x[-list_test,])
  pred_y <- as.matrix(y[-list_test,]); log(pred_y) -> pred_y
  
  print(rep)
  
  alpha=0.5
  tmp_results <- as.data.frame(matrix(ncol=5))
  colnames(tmp_results) <- c("rep", "alpha", "lambda", "r2_stand")
  
  tmp_results$rep <- rep; tmp_results$alpha <- alpha
  
  model<-cv.glmnet(test_x,test_y,nfolds=20,alpha=alpha,standardize=T)
  tmp_results$lambda <- model$lambda.min
  predicted <- predict(model, newx=pred_x, s="lambda.min")
  tmp_results$r2_stand <- summary(lm(predicted ~ pred_y))$r.squared
  weights<-t(as.matrix(unlist(coef(model,lambda="lambda.min"))[,1], nrow=1, ncol=30))
  
  rm(model, predicted)
  rbind(results_hist, tmp_results) -> results_hist; rm(tmp_results)
  rbind(betas_hist, weights) -> betas_hist; rm(weights)
  
  rm(alpha)
  rm(pred_x, pred_y, test_x, test_y)
  print("done")
}; rm(rep)

###############################
## clean up
###############################

rm(x, y, list_test, including_windows)

## create features (Supplementary Table)
as.data.frame(matrix(colnames(betas))) -> features
features$mean <- colMeans(betas)
features$historic <- colMeans(betas_hist)
features$recent <- colMeans(betas_recent)
features <- features[order(abs(features$mean), abs(features$recent), decreasing = T),]
write.table(features, "./glmnet_windows.txt", row.names=F, col.names=T, quote=F, sep = "\t")


###############################
## save output
###############################
save.image("./glmnet_model_results.RData")
b2 <- betas[, colSums(betas == 0)/nrow(betas) < 1]
b2[,-1] -> b2
beta_range <- as.data.frame(matrix(ncol=4, nrow=ncol(b2)))
colnames(beta_range) <- c('predictor', 'lower', 'upper', 'estimate')
to_plot <- melt(b2)

ggplot(data = to_plot, aes(x=Var2, y=value)) + 
  geom_boxplot() +
  coord_flip() +
  geom_violin() + 
  theme_set(theme_light())


