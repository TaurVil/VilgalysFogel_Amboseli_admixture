##### Initial data files #####
library(plyr); library(dplyr); library(data.table); library(parallel)
distance <- 2.5e5

# Start with human chromosome lengths (we used the fasta index file hg37) to get non-overlapping 250 kb windows of the genome
# Get matrix of Neanderthal ancestry for hg19 from Steinrucken et al. 

##### Get windows for human genome ######
lengths <- read.delim("~/genomes/hg37/hg37.fa.fai", header=F)

empty_windows <- NULL 
for (i in 1:22) {
  max <- lengths$V2[i]
  rows <- ceiling(max/distance); print(i); print(rows)
  tmp <- as.data.frame(matrix(ncol=3, nrow=rows))
  colnames(tmp) <- c("chr", "start", "end")
  tmp$chr <- lengths$V1[i]
  tmp$start <- (0:(rows-1))*distance
  tmp$end <- (1:rows)*distance
  tmp <- subset(tmp, tmp$start >= 5e4 & tmp$end <= (max-5e4))
  rbind(empty_windows,tmp) -> empty_windows; rm(tmp, max, rows)
}; rm(i)
rm(lengths)

##### Get Neanderthal ancestry ######
calc_mean_ancestry_per_window <- function(window) {
  subset(t2, t2$start <= walker$end[window] & t2$end > walker$start[window]) -> tmp2
  tmp2$start[tmp2$start < walker$start[window]] <- walker$start[window]
  tmp2$end[tmp2$end > walker$end[window]] <- walker$end[window]
  tmp2$l <- tmp2$end - tmp2$start
  return(sum(tmp2$l)/distance)
}

## list of all individuals published by Steinrucken et al.
individuals <- NULL
individuals$"CEU" = c('NA06984', 'NA06986', 'NA06989', 'NA06994', 'NA07000', 'NA07037', 'NA07048', 'NA07051', 'NA07056', 'NA07347', 'NA07357', 'NA10847', 'NA10851', 'NA11829', 'NA11830', 'NA11831', 'NA11843', 'NA11892', 'NA11893', 'NA11894', 'NA11919', 'NA11920', 'NA11930', 'NA11931', 'NA11932', 'NA11933', 'NA11992', 'NA11993', 'NA11994', 'NA11995', 'NA12003', 'NA12004', 'NA12006', 'NA12043', 'NA12044', 'NA12045', 'NA12046', 'NA12058', 'NA12144', 'NA12154', 'NA12155', 'NA12249', 'NA12272', 'NA12273', 'NA12275', 'NA12282', 'NA12283', 'NA12286', 'NA12287', 'NA12340', 'NA12341', 'NA12342', 'NA12347', 'NA12348', 'NA12383', 'NA12399', 'NA12400', 'NA12413', 'NA12489', 'NA12546', 'NA12716', 'NA12717', 'NA12718', 'NA12748', 'NA12749', 'NA12750', 'NA12751', 'NA12761', 'NA12763', 'NA12775', 'NA12777', 'NA12778', 'NA12812', 'NA12814', 'NA12815', 'NA12827', 'NA12829', 'NA12830', 'NA12842', 'NA12843', 'NA12872', 'NA12873', 'NA12874', 'NA12889', 'NA12890')
individuals$"CHB" = c('NA18525', 'NA18526', 'NA18527', 'NA18528', 'NA18530', 'NA18532', 'NA18534', 'NA18535', 'NA18536', 'NA18537', 'NA18538', 'NA18539', 'NA18541', 'NA18542', 'NA18543', 'NA18544', 'NA18545', 'NA18546', 'NA18547', 'NA18548', 'NA18549', 'NA18550', 'NA18552', 'NA18553', 'NA18555', 'NA18557', 'NA18558', 'NA18559', 'NA18560', 'NA18561', 'NA18562', 'NA18563', 'NA18564', 'NA18565', 'NA18566', 'NA18567', 'NA18570', 'NA18571', 'NA18572', 'NA18573', 'NA18574', 'NA18576', 'NA18577', 'NA18579', 'NA18582', 'NA18592', 'NA18593', 'NA18595', 'NA18596', 'NA18597', 'NA18599', 'NA18602', 'NA18603', 'NA18605', 'NA18606', 'NA18608', 'NA18609', 'NA18610', 'NA18611', 'NA18612', 'NA18613', 'NA18614', 'NA18615', 'NA18616', 'NA18617', 'NA18618', 'NA18619', 'NA18620', 'NA18621', 'NA18622', 'NA18623', 'NA18624', 'NA18626', 'NA18627', 'NA18628', 'NA18630', 'NA18631', 'NA18632', 'NA18633', 'NA18634', 'NA18635', 'NA18636', 'NA18637', 'NA18638', 'NA18639', 'NA18640', 'NA18641', 'NA18642', 'NA18643', 'NA18645', 'NA18647', 'NA18740', 'NA18745', 'NA18747', 'NA18748', 'NA18749', 'NA18757')
individuals$"CHS" = c('HG00403', 'HG00404', 'HG00406', 'HG00407', 'HG00418', 'HG00419', 'HG00421', 'HG00422', 'HG00427', 'HG00428', 'HG00436', 'HG00437', 'HG00442', 'HG00443', 'HG00445', 'HG00446', 'HG00448', 'HG00449', 'HG00451', 'HG00452', 'HG00457', 'HG00458', 'HG00463', 'HG00464', 'HG00472', 'HG00473', 'HG00475', 'HG00476', 'HG00478', 'HG00479', 'HG00500', 'HG00501', 'HG00512', 'HG00513', 'HG00524', 'HG00525', 'HG00530', 'HG00531', 'HG00533', 'HG00534', 'HG00536', 'HG00537', 'HG00542', 'HG00543', 'HG00556', 'HG00557', 'HG00559', 'HG00560', 'HG00565', 'HG00566', 'HG00577', 'HG00578', 'HG00580', 'HG00581', 'HG00583', 'HG00584', 'HG00589', 'HG00590', 'HG00592', 'HG00593', 'HG00595', 'HG00596', 'HG00607', 'HG00608', 'HG00610', 'HG00611', 'HG00613', 'HG00614', 'HG00619', 'HG00620', 'HG00625', 'HG00626', 'HG00628', 'HG00629', 'HG00634', 'HG00635', 'HG00650', 'HG00651', 'HG00653', 'HG00654', 'HG00656', 'HG00657', 'HG00662', 'HG00663', 'HG00671', 'HG00672', 'HG00683', 'HG00684', 'HG00689', 'HG00690', 'HG00692', 'HG00693', 'HG00698', 'HG00699', 'HG00701', 'HG00702', 'HG00704', 'HG00705', 'HG00707', 'HG00708')
n <- length(unique(individuals$CEU))

## Neanderthal ancestry data are split by chromosome, individual, and haplotype
## for each chromosome, get the relevant windows and 
## then for each individual load in data for that chromosome (split between two haplotypes)
## calculate the mean Neanderthal ancestry for each haplotype and window
Neand_ancestry <- NULL; 
for (chromosome in 1:22) { #chromosome=1
  c2 <- paste("chr",chromosome, sep="")
  win <- subset(empty_windows, empty_windows$chr == c2)
  max_i <- nrow(win)
  walker <- as.data.frame(matrix(ncol=sum(3,2*n), nrow=max_i) )
  win[,1:3] -> walker[,1:3]; colnames(win)[1:3] -> colnames(walker)[1:3]; rm(win)
  
  print(paste("chrom", chromosome)); print(max_i)
  
  z <- 4
  for (id in 1:n) { #id=1
    nom <- unique(individuals$CEU)[id]
    t2 <- read.delim(paste("~/my_genomes/hg37/Steinrucken2018_NeanderthalAncestry/CEU_lax_", c2, "/", nom, ".",chromosome,".0.bed", sep=""), header=F)
    colnames(t2) <- c('chr', 'start', 'end')
    walker[,z] <- do.call("c",mclapply(1:max_i, FUN = calc_mean_ancestry_per_window))
    colnames(walker)[z] <- paste(nom,"_0",sep="")
    z <- z+1
    t2 <- read.delim(paste("~/my_genomes/hg37/Steinrucken2018_NeanderthalAncestry/CEU_lax_", c2, "/", nom, ".",chromosome,".1.bed", sep=""), header=F)
    colnames(t2) <- c('chr', 'start', 'end')
    walker[,z] <- do.call("c",mclapply(1:max_i, FUN = calc_mean_ancestry_per_window))
    colnames(walker)[z] <- paste(nom,"_1",sep="")
    z <- z+1
    
    if (ceiling(id/45) == id/45) {print(paste(nom, "done"))}
  }; rm(id, max_i, t2, nom)
  rbind(Neand_ancestry,walker) -> Neand_ancestry; rm(walker)
}; rm(chromosome)
rm(c2, z, n, calc_mean_ancestry_per_window, individuals)
empty_windows -> all
all$mean_ancestry <- rowMeans(Neand_ancestry[,-c(1:3)])

##### Get genomic features ######

## Get B values ## these were hg18
## moved to hg19/37 using http://grch37.ensembl.org/Homo_sapiens/Tools/AssemblyConverter?db=core
## due to file sizes, we lifted them over in 3 parts (newB1-3)
calc_mean_B <- function(window) {
  subset(b, b$start <= tmp$end[window] & b$end > tmp$start[window]) -> tmp2
  tmp2$start[tmp2$start < tmp$start[window]] <- tmp$start[window]
  tmp2$end[tmp2$end > tmp$end[window]] <- tmp$end[window]
  tmp2$l <- tmp2$end - tmp2$start
  return(sum(tmp2$Bvalue*tmp2$l)/sum(tmp2$l))
}

all_b <- rbind(read.delim("~/my_genomes/hg37/newB1.bed", header=F), read.delim("~/my_genomes/hg37/newB2.bed", header=F), read.delim("~/my_genomes/hg37/newB3.bed", header=F))
colnames(all_b) <- c('chr', 'start', 'end', 'Bvalue')

all2 <- NULL; for (k in 1:22) {
  b <- subset(all_b, all_b$chr == k) 
  b$chr <- paste("chr",k,sep="")

  tmp <- subset(all, all$chr == paste("chr",k,sep=""))
  tmp$B <- do.call("c",mclapply(1:nrow(tmp), FUN = calc_mean_B))
  
  all2 <- rbind(all2, tmp); print(k) ; rm(tmp)
  rm(b)
}
all2 -> all; rm(all2, calc_mean_B, k, all_b, empty_windows)

## Get recombination (hg19, HapMap). rate is cM per Mb
calc_mean_rcr <- function(window) {
  end <- tmp$end[window]
  st <- tmp$start[window]
  k2 <- rcr22$start <= end & rcr22$end > st
  
  rcr22[c(k2==T),] -> tmp2
  tmp2$start[tmp2$start < tmp$start[window]] <- tmp$start[window]
  tmp2$end[tmp2$end > tmp$end[window]] <- tmp$end[window]
  tmp2$l <- tmp2$end - tmp2$start
  return(sum(tmp2$recombination*tmp2$l)/sum(tmp2$l))
}
all2 <- NULL; for (k in 1:22) {
  name=paste("~/my_genomes/hg37/HapMap_Recombination/genetic_map_GRCh37_chr",k,".txt",sep="")
  rcr22 <- fread(name)
  colnames(rcr22)[2] <- 'start'
  rcr22$end <- c(rcr22$start[-1], max(rcr22$start)+1)
  colnames(rcr22)[3] <- 'recombination'
  
  tmp <- subset(all, all$chr == paste("chr",k,sep=""))
  tmp$recombination <- do.call("c",mclapply(1:nrow(tmp), FUN = calc_mean_rcr))
  
  all2 <- rbind(all2, tmp); print(k) ; rm(tmp)
  rm(name, rcr22)
}; all2 -> all; rm(all2, k, calc_mean_rcr)

## Get derived archaic variants from http://ftp.eva.mpg.de/neandertal/ 
archaic_derived <- fread("~/my_genomes/hg37/ArchaicDerived_SNC_bothgq30.all_combined_maxsco_ranked.tsv")

calc_mean_fixed <- function(window) {
  end <- tmp$end[window]
  st <- tmp$start[window]
  k2 <- ad$Pos <= end & ad$Pos > st
  
  ad[c(k2==T),] -> tmp2
  return(nrow(tmp2))
}
all2 <- NULL; for (k in 1:22) {
  ad <- subset(archaic_derived, archaic_derived$`#Chrom` == k)
  tmp <- subset(all, all$chr == paste("chr",k,sep=""))
  tmp$fixed <- do.call("c",mclapply(1:nrow(tmp), FUN = calc_mean_fixed))
  
  all2 <- rbind(all2, tmp); print(k) ; rm(tmp, ad)
}; all2 -> all; rm(all2, k, calc_mean_fixed)
rm(archaic_derived)

## for each window, let's get the number of callable bases with Neanderthal genotypes
## To identify fixed differences, we will look for sites that don't vary in Neanderthals but are fixed in African populations of the 1000 genomes project
archaic_derived <- fread("~/my_genomes/hg37/Neand_and_human_frequencies/Neanderthal_fixed_variants_not_AFR.txt")
calc_mean_fixed <- function(window) {
  end <- tmp$end[window]
  st <- tmp$start[window]
  k2 <- ad$pos <= end & ad$pos > st
  
  ad[c(k2==T),] -> tmp2
  return(nrow(tmp2))
}
calc_possible <- function(window) {
  end <- tmp$end[window]
  st <- tmp$start[window]
  k2 <- testable$start <= end & testable$end > st
  
  testable[c(k2==T),] -> tmp2
  tmp2$start[tmp2$start <= st] <- st
  tmp2$end[tmp2$end >= end] <- end
  tmp2$l <- tmp2$end - tmp2$start
  return(sum(tmp2$l))
}

all2 <- NULL; for (k in 1:22) {
  ad <- subset(archaic_derived, archaic_derived$chr == k)
  testable <- fread(paste("~/my_genomes/hg37/Neand_and_human_frequencies/tmp.",k,".intervals",sep=""))
  testable$V2 <- testable$V2 + 1;colnames(testable) <- c('start','end')
  
  tmp <- subset(all, all$chr == paste("chr",k,sep=""))
  tmp$fixed_N3 <- do.call("c",mclapply(1:nrow(tmp), FUN = calc_mean_fixed))
  tmp$testable <- do.call("c",mclapply(1:nrow(tmp), FUN = calc_possible))
  
  
  all2 <- rbind(all2, tmp); print(k) ; rm(tmp, ad)
}; all2 -> all; rm(all2, k, calc_mean_fixed)
rm(testable, archaic_derived, calc_possible)

save.image("./windows.human.250kb.RData")
