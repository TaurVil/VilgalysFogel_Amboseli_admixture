## model expression data
library('EMMREML'); library('minpack.lm'); library(gap); library(cobs)

dataset="eLife"
load(paste("./expression_data_for_models.", dataset, ".RData", sep=""))

#######################################
## run models using EMMREML
#######################################
subset(residuals, row.names(residuals) %in% row.names(loc_anc)) -> data; dim(data)
data[order(row.names(data)),] -> data; sum(!(row.names(data) == row.names(loc_anc)))
output1 <- as.data.frame(matrix(ncol=7, nrow=nrow(data))); colnames(output1) <- c('gene', 'beta_local', 'p_local', 'beta_global', 'p_global', 'se_local', 'se_global')
rownames(data) -> output1[,1]
for (i in c(1:nrow(data))) {
  t(as.matrix(loc_anc[i,])) -> temp_ancestry
  t(as.matrix(data[i,])) -> temp_ge
  design <- cbind(rep(1,nrow(temp_ancestry)), temp_ancestry, metadata$mean_ancestry, metadata$age, metadata$sex)
  
  emma=emmreml(y=as.vector(temp_ge),X=design,Z=as.matrix(Z),K=as.matrix(gt_cov2), varbetahat=T, varuhat=T, PEVuhat=T,test=T)
  output1[i,2] <- emma$betahat[2]; output1[i,3]=emma$pvalbeta[2,8] # local ancestry
  output1[i,4] <- emma$betahat[3]; output1[i,5]=emma$pvalbeta[3,8] # mean ancestry
  output1[i,6:7] <- emma$varbetahat[2:3] # SE of local and mean ancestry
  rm(temp_ancestry, temp_ge, design)
}; rm(i, emma)
output1 -> results
# hist(results$p_local, breaks=200)

#######################################
## permute local ancestry
#######################################
set.seed(42)
output1 <- as.data.frame(matrix(ncol=7, nrow=nrow(data))); colnames(output1) <- c('gene', 'beta_local', 'p_local', 'beta_global', 'p_global', 'se_local', 'se_global')
rownames(data) -> output1[,1]
for (i in c(1:nrow(data))) {
  t(as.matrix(sample(loc_anc[i,]))) -> temp_ancestry
  t(as.matrix(data[i,])) -> temp_ge
  design <- cbind(rep(1,nrow(temp_ancestry)), temp_ancestry, metadata$mean_ancestry, metadata$age, metadata$sex)
  
  emma=emmreml(y=as.vector(temp_ge),X=design,Z=as.matrix(Z),K=as.matrix(gt_cov2), varbetahat=T, varuhat=T, PEVuhat=T,test=T)
  output1[i,2] <- emma$betahat[2]; output1[i,3]=emma$pvalbeta[2,8]
  output1[i,4] <- emma$betahat[3]; output1[i,5]=emma$pvalbeta[3,8]
  output1[i,6:7] <- emma$varbetahat[2:3] # SE in the first two betas
  rm(temp_ancestry, temp_ge, design)
}; output1 -> perm_results
for (k in 1:49) {
  output1 <- as.data.frame(matrix(ncol=7, nrow=nrow(data))); colnames(output1) <- c('gene', 'beta_local', 'p_local', 'beta_global', 'p_global', 'se_local', 'se_global')
  rownames(data) -> output1[,1]
  for (i in c(1:nrow(data))) {
    t(as.matrix(sample(loc_anc[i,]))) -> temp_ancestry
    t(as.matrix(data[i,])) -> temp_ge
    design <- cbind(rep(1,nrow(temp_ancestry)), temp_ancestry, metadata$mean_ancestry, metadata$age, metadata$sex)
    
    emma=emmreml(y=as.vector(temp_ge),X=design,Z=as.matrix(Z),K=as.matrix(gt_cov2), varbetahat=T, varuhat=T, PEVuhat=T,test=T)
    output1[i,2] <- emma$betahat[2]; output1[i,3]=emma$pvalbeta[2,8]
    output1[i,4] <- emma$betahat[3]; output1[i,5]=emma$pvalbeta[3,8]
    output1[i,6:7] <- emma$varbetahat[2:3] # SE in the first two betas
    rm(temp_ancestry, temp_ge, design)
  }
  rbind(output1, perm_results) -> perm_results
  print(nrow(perm_results)/nrow(output1))
}; rm(emma, k)

#######################################
## permute global ancestry
#######################################
set.seed(42)
output1 <- as.data.frame(matrix(ncol=7, nrow=nrow(data))); colnames(output1) <- c('gene', 'beta_local', 'p_local', 'beta_global', 'p_global', 'se_local', 'se_global')
rownames(data) -> output1[,1]
for (i in c(1:nrow(data))) {
  t(as.matrix((loc_anc[i,]))) -> temp_ancestry
  t(as.matrix(data[i,])) -> temp_ge
  design <- cbind(rep(1,nrow(temp_ancestry)), temp_ancestry, sample(metadata$mean_ancestry), metadata$age, metadata$sex)
  
  emma=emmreml(y=as.vector(temp_ge),X=design,Z=as.matrix(Z),K=as.matrix(gt_cov2), varbetahat=T, varuhat=T, PEVuhat=T,test=T)
  output1[i,2] <- emma$betahat[2]; output1[i,3]=emma$pvalbeta[2,8]
  output1[i,4] <- emma$betahat[3]; output1[i,5]=emma$pvalbeta[3,8]
  output1[i,6:7] <- emma$varbetahat[2:3] # SE in the first two betas
  rm(temp_ancestry, temp_ge, design)
}; output1 -> perm_global
for (k in 1:49) {
  output1 <- as.data.frame(matrix(ncol=7, nrow=nrow(data))); colnames(output1) <- c('gene', 'beta_local', 'p_local', 'beta_global', 'p_global', 'se_local', 'se_global')
  rownames(data) -> output1[,1]
  for (i in c(1:nrow(data))) {
    t(as.matrix((loc_anc[i,]))) -> temp_ancestry
    t(as.matrix(data[i,])) -> temp_ge
    design <- cbind(rep(1,nrow(temp_ancestry)), temp_ancestry, sample(metadata$mean_ancestry), metadata$age, metadata$sex)
    
    emma=emmreml(y=as.vector(temp_ge),X=design,Z=as.matrix(Z),K=as.matrix(gt_cov2), varbetahat=T, varuhat=T, PEVuhat=T,test=T)
    output1[i,2] <- emma$betahat[2]; output1[i,3]=emma$pvalbeta[2,8]
    output1[i,4] <- emma$betahat[3]; output1[i,5]=emma$pvalbeta[3,8]
    output1[i,6:7] <- emma$varbetahat[2:3] # SE in the first two betas
    rm(temp_ancestry, temp_ge, design)
  }
  rbind(output1, perm_global) -> perm_global
  print(nrow(perm_global)/nrow(output1))
}; rm(emma, k)

#######################################
## permutation-based FDR correction, using Joaquin's code from Snyder-Mackler et al. 2016
#######################################

library(cobs)
perm.fdr=function(input_df,perm_df,Pvals_col_name,name){
  
  pvals_index=which(colnames(input_df)==Pvals_col_name)
  ro<-input_df[order(input_df[,pvals_index]),]
  p_obs <- data.frame(pvalue=ro[,pvals_index])
  p_vector<-matrix(as.matrix(perm_df),ncol=1)
  p_vector=data.frame(p_vector[order(p_vector)])
  
  F<-p_obs[,1]
  F_o<-p_obs[,1]
  pi_hat<-p_obs[,1]
  
  j=1
  observed<-length(p_obs[,1])
  randoms<-length(p_vector[,1])
  
  for(i in 1:observed)
  {
    repeat
    {
      if((p_vector[j,1]<p_obs[i,1])&j<randoms){j<-j+1}else{break}
    }
    F[i]=i/observed
    F_o[i]=(j-1)/randoms
    if(F_o[i]<1){pi_hat[i]=(1-F[i])/(1-F_o[i])}else{pi_hat[i]=1}
  }
  tabla <-data.frame(pi_hat,pval=p_obs[,1])
  
  tabla[1,]=c(1,0)
  last_percentile_average=mean(tabla$pi_hat[as.integer(min((length(tabla[,1])*0.99),(nrow(tabla)-1)):length(tabla[,1]))])
  tabla[nrow(tabla),]=c(last_percentile_average,1)
  constraint_matrix=as.matrix(data.frame(c(0,2),c(0,1),c(1,0)))
  f_hat<-suppressWarnings(cobs(tabla$pval,tabla$pi_hat,constraint="convex",pointwise=constraint_matrix,maxiter=1000,print.warn=FALSE,print.mesg=FALSE))
  
  f_hat_serie=f_hat$fitted
  pi_o=f_hat_serie[length(f_hat_serie)]
  pi_o=min(pi_o,1)
  pi_o=max(pi_o,0)
  
  Fdr_ST_perm=pi_o*F_o/F
  
  for(i in 1:length(p_obs[,1]))
  {
    Fdr_ST_perm[i]=pi_o*F_o[i]/F[i]
    if(i>1)
    {
      for(j in 1:(i-1))
      {
        if(Fdr_ST_perm[i-j]>Fdr_ST_perm[i]){Fdr_ST_perm[i-j]=Fdr_ST_perm[i]}else{break}
      }
    }
    if(Fdr_ST_perm[i]>1)  Fdr_ST_perm[i]=1
  }
  
  fdrs_df <-data.frame(ro,q_ST_perm=Fdr_ST_perm)
  rownames(fdrs_df)=rownames(ro)
  colnames(fdrs_df)[ncol(fdrs_df)]=paste0("fdr_",name)
  
  return(fdrs_df)
}

#######################################
## clean and save
#######################################
rm(i, lane, Z, residuals, loc_anc, output1, gt_cov2, data, metadata)

if (dataset == "TruCulture") {
   gene_lengths -> genes_TC; rm(gene_lengths)
   results -> TC; rm(results)
   perm_results -> perm_TC; rm(perm_results)
   perm_global -> perm_global_TC; rm(perm_global)
}
if (dataset == "eLife") {
   gene_lengths -> genes_eLife; rm(gene_lengths)
   results -> eLife; rm(results)
   perm_results -> perm_eLife; rm(perm_results)
   perm_global -> perm_global_eLife; rm(perm_global)
}

save.image(paste("./expression_results.", dataset, ".RData", sep=""))
 
