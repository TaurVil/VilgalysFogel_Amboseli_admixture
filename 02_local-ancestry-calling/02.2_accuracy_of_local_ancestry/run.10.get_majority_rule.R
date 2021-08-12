#!/usr/bin/env Rscript
#SBATCH --get-user-env

library(data.table)
print(version$version.string)
print(paste("data.table version", packageVersion("data.table")))

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

#get data for the same
fread("./raw_calls/DISTANCE.DIFFERENCE.INDIV.COVERAGE.txt") -> indv
colnames(indv) <- c('chrom', 'snp', 'call')
indv$indiv <- "INDIV"
read.delim("~/genomes/panubis1/Panubis1.0.fa.fai", header=F) -> chroms; subset(chroms, chroms$V1 %in% unique(indv$chrom)) -> chroms

#We're going to get ancestry tracts for each chrom individually (going through the intermediary of majority rule)
#Then, at the end of each chromosome, we'll merge the tracts for that chromosome with the building list of all tracts for that individual

# Value sets the window size
# Mode sets the minimum percentage for a call to be made using majority rule
# min_n sets the minimum number of SNPs within the window to make a call
# exclude refers to how much of the start/end of each chromosome to ignore
value=DISTANCE
mode_n=0.5
min_n=20
exclude=50000

for (j in 1:nrow(chroms)) {
  subset(indv, indv$chrom == chroms[j,1]) -> tmp
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
    ##remove sites where there is not a concensus call by majority rule, or not enough sites to make a confident call
    subset(tmp, i1$perc >= mode_n & i1$n >= min_n) -> tmp; subset(i1, i1$perc >= mode_n & i1$n >= min_n) -> i1
    if (nrow(i1) > 0) {
      ##Merge majority rule calls for this chromosome to the growing file.
      if (j==1) {cbind(tmp,i1)[tmp$snp >= exclude & tmp$snp <= chroms[j,2]-exclude,-c(5:6)] -> maj_rule} else {cbind(tmp,i1) -> te; rbind(maj_rule, te[tmp$snp >= exclude & tmp$snp <= chroms[j,2]-exclude,-c(5:6)]) -> maj_rule}
      
      ##Now let's move onto calling ancestry tracts now that we have the majority rules.
      cbind(tmp, i1) -> modes; rm(tmp); rm(i1); if (j > 1) {rm(te)}
      
      modes$nxt <- c(as.numeric(modes$snp[-1]),max(modes$snp))
      modes$nxt_chrom <- c(as.character(modes$chrom[-1]),modes$chrom[nrow(modes)])
      modes$nxt_state <- c(modes$mode[-1],modes$call[nrow(modes)])
      
      #Keep only sites where this site is different than the next one.
      ##aka, first row is now the first break point
      subset(modes, !(modes$mode == modes$nxt_state & modes$chrom == modes$nxt_chrom)) -> blocks
      blocks[,.(chrom, nxt_chrom, snp, nxt, mode, nxt_state)] -> blocks
      
      ##nxt_state is the one that defines that block
      ###Let's add the first block from position 1 to the first break pt
      first <- as.data.frame(t(matrix(c(blocks$chrom[1], blocks$chrom[1], 1, blocks$snp[1], modes$mode[1], blocks$mode[1]))))
      colnames(first) <- colnames(blocks); rbind(first, blocks) -> blocks; rm(first)
      
      last <- as.data.frame(matrix(c(as.character(chroms[j,1]), as.character(chroms[j,1]), max(modes$snp), chroms[j,2], modes$nxt_state[nrow(modes)], modes$nxt_state[nrow(modes)]), ncol=6)); as.numeric(as.character(last$V3)) -> last$V3; as.numeric(as.character(last$V4)) -> last$V4; as.numeric(as.character(last$V5)) -> last$V5;as.numeric(as.character(last$V6)) -> last$V6;
      colnames(last) <- colnames(blocks); rbind(blocks, last) -> blocks; rm(last)
      
      blocks$chrom <- blocks$nxt_chrom <- as.character(chroms[j,1])
      
      write.table(blocks, "./tmp2.DISTANCE.DIFFERENCE.INDIV.COVERAGE.txt", row.names=F, col.names=T, sep="\t")
      read.delim("./tmp2.DISTANCE.DIFFERENCE.INDIV.COVERAGE.txt") -> blocks
      
      blocks$brk <- (as.numeric(blocks$nxt) + as.numeric(blocks$snp))/2
      blocks$brk[1] <- 1
      
      ## We want the first break to be at the first snp
      chroms[j,2] -> blocks$brk[nrow(blocks)]
      ## The last break is already set at the last snp
      
      blocks$nxt_brk <- c(as.numeric(blocks$brk[-1]),'end')
      #blocks$nxt_brk[nrow(blocks)-1] <- blocks$snp[nrow(blocks)]
      blocks[-nrow(blocks),] -> blocks
      
      blocks$length <- as.numeric(blocks$nxt_brk) - as.numeric(blocks$brk)
      #blocks$length[1] <- NA; blocks$length[nrow(blocks)-1] <- NA
      
      ## Removed uncertainty because we really don't need this now that I'm not doing as much methods testing.
      #uncertainty of the break points to each side
      blocks$u_prev <- as.numeric(blocks$nxt) - as.numeric(blocks$snp)
      blocks$u_prev[1] <- NA
      blocks$u_next <- c(as.numeric(blocks$u_prev[-1]),NA)
      
      #Blocks starts at 1 and goes to the end of the chromosome.
      
      t <- blocks[, c(1, 7, 8, 6, 9:11)]
      colnames(t) <- c("chrom", "start", "end", "state", "length", "prev_inferred", "after_inferred")
      t$name <- "INDIV"
      #t$length2 <- as.numeric(t$end) - as.numeric(t$start)
      t$end <- as.numeric(as.character(t$end))
      t$start <- as.numeric(as.character(t$start))
      
      # we'll fix the extrmes of t so we crop the first and last 'exclude'bp from the chromosome
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

write.table(tracts, "./tracts/lclae.DISTANCE.DIFFERENCE.majruleDISTANCE.min50.n20.exclude50kb.INDIV.COVERAGE.tracts.txt", row.names=F, col.names=T, quote=F, sep="\t")
