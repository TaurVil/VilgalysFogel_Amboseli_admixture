#######################################
## load in results. 
## define function for FDR calculation
#######################################
dataset="eLife"; load(paste("./RESULTS/expression_piecewise.", dataset, ".RData", sep=""))
dataset="TruCulture"; load(paste("./RESULTS/expression_piecewise.", dataset, ".RData", sep=""))
rm(dataset)

library(cobs); perm.fdr=function(input_df,perm_df,Pvals_col_name,name){
  
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
## filter for sites with enough coverage and calculate FDR
#######################################

res_eLife <- res_eLife[apply(cbind(res_eLife$anu, res_eLife$yel, res_eLife$het), 1, min)>=10,]
res_TC <- res_TC[apply(cbind(res_TC$anu, res_TC$yel, res_TC$het), 1, min)>=10,]
perm_eLife <- perm_eLife[apply(cbind(perm_eLife$anu, perm_eLife$yel, perm_eLife$het), 1, min)>=10,]
perm_TC <- perm_TC[apply(cbind(perm_TC$anu, perm_TC$yel, perm_TC$het), 1, min)>=10,]

res_eLife <- res_eLife[!is.na(res_eLife$slope2p),]
perm_eLife <- perm_eLife[!is.na(perm_eLife$slope2p),]
res_TC <- res_TC[!is.na(res_TC$slope2p),]
perm_TC <- perm_TC[!is.na(perm_TC$slope2p),]
res_TC <- perm.fdr(res_TC, perm_TC$slope2p, "slope2p", "10perm")
res_eLife <- perm.fdr(res_eLife, perm_eLife$slope2p, "slope2p", "10perm")

sum(res_eLife$fdr_10perm < 0.1)
sum(res_TC$fdr_10perm < 0.1)
res_TC[1:sum(res_TC$fdr_10perm < 0.1),]

#######################################
## plot qq plots
## only a weak enrichment in the larger dataset
#######################################

par(mfrow=c(1,2)) 
t1 <- -log10(res_TC$slope2p); t2 <- -log10(perm_TC$slope2p)
qqplot(t2, t1, xlab="permuted", ylab="observed", main="PBMCs\n(n=94 samples)"); abline(a=0,b=1,col='red')
t1 <- -log10(res_eLife$slope2p); t2 <- -log10(perm_eLife$slope2p)
qqplot(t2, t1, xlab="permuted", ylab="observed", main="whole blood\n(n=63 samples)"); abline(a=0,b=1,col='red')

