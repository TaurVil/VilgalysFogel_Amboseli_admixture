#!/usr/bin/env Rscript
#SBATCH --get-user-env

library(data.table)

print(version$version.string)
print(packageVersion("data.table"))

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

#get data for the sample
fread("./INDIV.35kb.d2.masked.SNPRCref.txt") -> indv
colnames(indv) <- c('chrom', 'snp', 'call', 'indiv')

# Create list of chromosomes to run
# updated to Panubis1.0, and only return chromosomes for which LCLAE was run
read.delim("./DATA/Panubis1.0.fa.fai", header=F) -> chroms
#subset(chroms, chroms$V1 %in% unique(indv$chrom)) -> chroms # This needs to be commented out when dealing with the "chr" in the chromosome names, because otherwise they don't match at all

# We're going to get ancestry tracts for each chromosome individually (going through the intermediary of majority rule) 
# Then, at the end of each chromosome, we'll merge the tracts for that chromosome with the building list of all tracts for that individual

# Value sets the window size
# Mode sets the minimum proportion of calls for a call to be made using majority rule
value=35000 #window size for majority rules
mode_n=0.5 #minimum proportion of calls for majority rule
min_n=20 #minimum number of nearby ancestry-informative sites to retain an ancestry call
exclude=50000 ## new option to not return the first/last X bp of each chromosome

##Run a for-loop for each chromosome
for (j in 1:20) {  ## Chromosomes 1-20
  subset(indv, indv$chrom == j) -> tmp
  ## -1 calls occur when LCLAE is not able to assign a most likely ancestry for a ancestry-informative marker (e.g. two ancestry states are equally likely)
  ## Alternatively, these sites can be removed before using this script using `grep -v '-1' ./INDIV.35kb.d2.masked.SNPRCref.txt` {tpv: this is the appraoch I've used recently}
  subset(tmp, !(tmp$call == -1)) -> tmp
  ## If we don't make any calls, just skip to the next chromosome
if (nrow(tmp) > 0) {
  ## Index and setkeys for match. Then run match to get sets of sites within $VALUE of each SNP. Then get the mode of all calls, the total number of calls, and the number of calls which were the mode. 
  tmp[,loc_Dummy := snp]; tmp[,.(snp, call)] -> tmp2
  tmp2[,loc_Plus100 := snp + value]; tmp2[,loc_Minus100 := snp - value]
  setkey(tmp,snp,loc_Dummy); setkey(tmp2,loc_Minus100, loc_Plus100)
  print(nrow(tmp))
  print(paste("Now doing matches for", j, "....", sep=""))
  Matches <- foverlaps(tmp[,.(snp, loc_Dummy)], tmp2[,.(loc_Minus100,loc_Plus100,call)])
  Matches[,.(n = .N, mode = getmode(call), n_mode=sum(call==getmode(call))), by = .(snp)] -> i1
  rm(Matches); gc(); rm(tmp2)
  print("done with matches!")
  i1$perc <- i1$n_mode/i1$n
  ## Remove sites where there is not a concensus call by majority rule (meating mode_n) or enough nearby ancestry informative sites (min_n)
  subset(tmp, i1$perc >= mode_n & i1$n >= min_n) -> tmp; subset(i1, i1$perc >= mode_n & i1$n >= min_n) -> i1
if (nrow(i1) > 0) {
  ## Merge majority rule calls for this chromosome to the growing file, removing the first and last X kb (set by "exclude"). 
  if (j==1) {cbind(tmp,i1)[tmp$snp >= exclude & tmp$snp <= chroms[j,2]-exclude,-c(5:6)] -> maj_rule} else {cbind(tmp,i1) -> te; rbind(maj_rule, te[tmp$snp >= exclude & tmp$snp <= chroms[j,2]-exclude,-c(5:6)]) -> maj_rule}
  
  ## Now let's move onto calling ancestry tracts now that we have the majority rule calls 
  cbind(tmp, i1) -> modes; rm(tmp); rm(i1); if (j > 1) {rm(te)}
  
  modes$nxt <- c(as.numeric(modes$snp[-1]),max(modes$snp))
  modes$nxt_chrom <- c(as.character(modes$chrom[-1]),modes$chrom[nrow(modes)])
  modes$nxt_state <- c(modes$mode[-1],modes$call[nrow(modes)])
  
  ## Keep only sites where this site is different than the next one. 
  ## AKA, the first row is now the first break point
  subset(modes, !(modes$mode == modes$nxt_state & modes$chrom == modes$nxt_chrom)) -> blocks
  blocks[,.(chrom, nxt_chrom, snp, nxt, mode, nxt_state)] -> blocks
  
  ## nxt_state is the one that defines that block
  ### Let's add the first block from position 1 to the first break pt
  first <- as.data.frame(t(matrix(c(blocks$chrom[1], blocks$chrom[1], 1, blocks$snp[1], modes$mode[1], blocks$mode[1]))))
  colnames(first) <- colnames(blocks); rbind(first, blocks) -> blocks; rm(first)

  last <- as.data.frame(matrix(c(as.character(chroms[j,1]), as.character(chroms[j,1]), max(modes$snp), chroms[j,2], modes$nxt_state[nrow(modes)], modes$nxt_state[nrow(modes)]), ncol=6)); as.numeric(as.character(last$V3)) -> last$V3; as.numeric(as.character(last$V4)) -> last$V4; as.numeric(as.character(last$V5)) -> last$V5
;as.numeric(as.character(last$V6)) -> last$V6;
  colnames(last) <- colnames(blocks); rbind(blocks, last) -> blocks; rm(last)

  blocks$chrom <- blocks$nxt_chrom <- as.character(chroms$V1[j]) # If you want the numbers, this should be `j` instead of `chroms$V1[j]`

  write.table(blocks, "./tmp2.INDIV.txt", row.names=F, col.names=T, sep="\t")
  read.delim("./tmp2.INDIV.txt") -> blocks
  
  blocks$brk <- (as.numeric(blocks$nxt) + as.numeric(blocks$snp))/2
  blocks$brk[1] <- 1  # We start at the first bp of the chromosome
  modes$snp[1] -> blocks$brk[1] # Which we'll assume has the same ancestry state as the first called SNP
 
  chroms[j,2]  -> blocks$brk[nrow(blocks)] # The last block should end at the end of the chromosome
  
  
  blocks$nxt_brk <- c(as.numeric(blocks$brk[-1]),'end')
  blocks[-nrow(blocks),] -> blocks  # last block no longer matters, we've already extended from the last SNP to the end of the chromosome

  blocks$length <- as.numeric(blocks$nxt_brk) - as.numeric(blocks$brk)
  #blocks$length[1] <- NA; blocks$length[nrow(blocks)-1] <- NA
  
  ## Removed uncertainty because we really don't need this now that I'm not doing as much methods testing. 
  #uncertainty of the break points to each side
  blocks$u_prev <- (as.numeric(blocks$nxt) - as.numeric(blocks$snp))/2
  blocks$u_prev[1] <- NA
  blocks$u_next <- c(as.numeric(blocks$u_prev[-1]),NA)
  ## blocks u_prev is the number of bases at the start of that tract which were inferred (i.e. before the first AIM with that ancestry call)
  ## blocks u_next is the same thing for the end of that tract
  
  #Blocks starts at the first SNP and goes to the last one. 
  
  t <- blocks[, c(1, 7, 8, 6, 9:11)]
  colnames(t) <- c("chrom", "start", "end", "state", "length", "prev_inferred", "after_inferred")
  t$name <- "INDIV"  
  #t$length2 <- as.numeric(t$end) - as.numeric(t$start)
  t$end <- as.numeric(as.character(t$end))
  t$start <- as.numeric(as.character(t$start))

  # We'll fix the extremes of t so we crop the first and last 'exclude'bp from the chromosome
  t <- subset(t, t$start <= (chroms[j,2]-exclude) & t$end >= exclude)
  t$start[t$start < exclude] <- exclude
  t$end[t$end > (chroms[j,2]-exclude)] <- (chroms[j,2]-exclude)
  t$length[c(1,nrow(t))] <- NA

  if (j==1) {t -> tracts} else {rbind(tracts,t) -> tracts}
  rm(t); rm(modes) 
  rm(blocks)
}
}
}

write.table(maj_rule, "INDIV_MajorityRule.d2.min50.masked.SNPRCref.txt", row.names=F, col.names=T, quote=F, sep="\t")
write.table(tracts, "INDIV_tracts.d2.35kb.masked.SNPRCref.txt", row.names=F, col.names=T, quote=F, sep="\t")
