## model non-additive effects in the expression data
library('EMMREML'); library('minpack.lm'); library(gap); library(cobs); library(segmented)

dataset="TruCulture"
load(paste("./DATA/expression_data_for_models.", dataset, ".RData", sep=""))

#######################################
## regress out the mixed effect, producing "resids" for the residuals
#######################################
subset(residuals, row.names(residuals) %in% row.names(loc_anc)) -> data; dim(data)
data[order(row.names(data)),] -> data; sum(!(row.names(data) == row.names(loc_anc)))
resids  <- data; resids[] <- NA
for (i in c(1:nrow(data))) {
  t(as.matrix(loc_anc[i,])) -> temp_ancestry
  t(as.matrix(data[i,])) -> temp_ge
  design <- cbind(rep(1,nrow(temp_ancestry)), temp_ancestry, metadata$mean_ancestry, metadata$age, metadata$sex)
  
  emma=emmreml(y=as.vector(temp_ge),X=design,Z=as.matrix(Z),K=as.matrix(gt_cov2), varbetahat=T, varuhat=T, PEVuhat=T,test=T)
  temp_ge - t(emma$uhat) -> resids[i,]
  
  rm(temp_ancestry, temp_ge, design)
}; rm(i, emma, Z)
rm(lane, gt_cov2)

#######################################
## fit a stepwise function
## Coerce individuals into nearest state (homozygous or heterozygous) to count the number of samples per state
## Consider only genes with at least 5 individuals from each ancestry state
#######################################

# Prepare variables
metadata$mean_ancestry -> Global.Ancestry; metadata$age -> Age; metadata$sex -> Sex

## for each site, run segmented model if there's sufficient ancestry information
res <- as.data.frame(matrix(ncol=13,nrow=nrow(resids))); row.names(res) <- row.names(resids)
colnames(res) <- c('logLik1', 'logLik2', 'slope1', 'slope2', 'slope1p', 'slope2p', "SE1", "SE2", "daviesP", "daviesStat", "anu", "yel", "het")
for (i in 1:nrow(res)) {
  t(as.matrix(loc_anc[i,])) -> temp_ancestry
  (as.matrix(resids[i,])) -> temp_ge
  
  sum(temp_ancestry >= 1.5) -> a -> res[i,11]; sum(temp_ancestry <= 0.5) -> y -> res[i,12]; sum(temp_ancestry > 0.5 & temp_ancestry < 1.5) -> het -> res[i,13]
  
  if (min(a, y, het) >= 5) {
    ## Simple linear model, 6 degrees freedom
    lin.mod <- lm((temp_ge)~(temp_ancestry) + Global.Ancestry + Age + Sex)
    logLik(lin.mod)[1] -> res[i,1]
    ## Run a segmented linear model
    segmented.mod <- segmented(lin.mod, seg.Z = ~temp_ancestry, psi=1, control=seg.control(it.max=0))
    summary(segmented.mod)$coefficients[2,1] -> res[i,3]
    summary(segmented.mod)$coefficients[6,1] -> res[i,4]
    summary(segmented.mod)$coefficients[2,4] -> res[i,5]
    summary(segmented.mod)$coefficients[6,4] -> res[i,6]
    summary(segmented.mod)$coefficients[2,2] -> res[i,7]
    summary(segmented.mod)$coefficients[6,2] -> res[i,8]
    logLik(segmented.mod)[1] -> res[i,2]
    d.test <- davies.test(lin.mod, seg.Z = ~temp_ancestry, k=100)
    res[i,9] <- d.test$p.value
    res[i,10] <- d.test$statistic[1]
    rm(d.test, lin.mod, segmented.mod)
    }
  
  rm(temp_ge, temp_ancestry, a, y, het)
}
res -> res_nonadditive

#######################################
## permutations of local ancestry in the stepwise function
#######################################
for (k in 1:10) {
  res <- as.data.frame(matrix(ncol=13,nrow=nrow(resids))); row.names(res) <- row.names(resids)
  colnames(res) <- c('logLik1', 'logLik2', 'slope1', 'slope2', 'slope1p', 'slope2p', "SE1", "SE2", "daviesP", "daviesStat", "anu", "yel", "het")
  for (i in 1:nrow(resids)) {
    t(as.matrix(sample(loc_anc[i,]))) -> temp_ancestry
    (as.matrix(resids[i,])) -> temp_ge
    
    sum(temp_ancestry >= 1.5) -> a -> res[i,11]; sum(temp_ancestry <= 0.5) -> y -> res[i,12]; sum(temp_ancestry > 0.5 & temp_ancestry < 1.5) -> het -> res[i,13]
    
    if (min(a, y, het) >= 10) {
      lin.mod <- lm((temp_ge)~(temp_ancestry) + Global.Ancestry + Age + Sex)
      logLik(lin.mod)[1] -> res[i,1]
      #tryCatch({
      segmented.mod <- segmented(lin.mod, seg.Z = ~temp_ancestry, psi=1, control=seg.control(it.max=0))
      summary(segmented.mod)$coefficients[2,1] -> res[i,3]
      summary(segmented.mod)$coefficients[6,1] -> res[i,4]
      summary(segmented.mod)$coefficients[2,4] -> res[i,5]
      summary(segmented.mod)$coefficients[6,4] -> res[i,6]
      summary(segmented.mod)$coefficients[2,2] -> res[i,7]
      summary(segmented.mod)$coefficients[6,2] -> res[i,8]
      logLik(segmented.mod)[1] -> res[i,2]
      #}, error=function(e){''})
      d.test <- davies.test(lin.mod, seg.Z = ~temp_ancestry, k=100)
      res[i,9] <- d.test$p.value
      res[i,10] <- d.test$statistic[1]
      rm(temp_ge, temp_ancestry, a, y, het)
    }
  }
  if (k==1) {res -> perm_nonadditive} else {rbind(res, perm_nonadditive) -> perm_nonadditive}
  rm(res, i); print(k)
}; rm(k)
rm(segmented.mod, lin.mod, d.test, Age, Global.Ancestry, Sex, ids, data, gene_lengths, residuals)

#######################################
## clean and save data file
#######################################
if (dataset == "TruCulture") {
  perm_nonadditive -> perm_TC; rm(perm_nonadditive)
  res_nonadditive -> res_TC; rm(res_nonadditive)
  resids -> data_TC; rm(resids)
  metadata -> meta_TC; rm(metadata)
  loc_anc -> la_TC; rm(loc_anc)
}
if (dataset == "eLife") {
  perm_nonadditive -> perm_eLife; rm(perm_nonadditive)
  res_nonadditive -> res_eLife; rm(res_nonadditive)
  resids -> data_eLife; rm(resids)
  metadata -> meta_eLife; rm(metadata)
  loc_anc -> la_eLife; rm(loc_anc)
}

save.image(paste("./RESULTS/expression_piecewise.", dataset, ".RData", sep=""))


