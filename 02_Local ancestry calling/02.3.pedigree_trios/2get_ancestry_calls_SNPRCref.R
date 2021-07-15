#!/usr/bin/env Rscript

load("local_ancestry_pedigree_trios_maskedSNPRCref.Rd", verbose=T)

# Run R script separately for each chromosome
tmp <- subset(positions, chrom=="CHROMOSOME") # only positions for focal chromosome
tmp2 <- subset(tracts2, chrom=="CHROMOSOME") # only tracts from focal chromosome

ancestry_calls <- as.data.frame(matrix(ncol=(nrow(indiv_list)+2), nrow=nrow(tmp)))
colnames(ancestry_calls)[1:2] <- c("chrom", "pos")
colnames(ancestry_calls)[3:(nrow(indiv_list)+2)] <- as.character(indiv_list[,1])
ancestry_calls$chrom <- tmp[,1]
ancestry_calls$pos <- tmp[,2]

for (k in 1:nrow(indiv_list)) {
    
    tmp3 <- tmp2[tmp2$indiv==indiv_list[k,1],]
    
  for (j in 1:nrow(tmp)) {
  
  tmp4 <- subset(tmp3, tmp3$start <= ancestry_calls$pos[j] & tmp3$end >= ancestry_calls$pos[j])
  
      if (nrow(tmp4) == 1) { ancestry_calls[j,k+2] <- tmp4$state }
      if (nrow(tmp4) > 1) { ancestry_calls[j,k+2] <- "boundary"}
      if (nrow(tmp4) == 0) { ancestry_calls[j,k+2] <- "not_in_tract"}
  
  }
 print(k)
     
}


write.table(ancestry_calls,file="ancestry_calls_maskedSNPRCref_pedigree_trios_35kbpos_CHROMOSOME.txt",quote=F,row.names=F)
print("done")
q(save="no")
