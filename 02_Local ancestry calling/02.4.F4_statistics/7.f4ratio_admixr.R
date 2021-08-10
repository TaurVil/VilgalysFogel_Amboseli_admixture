#!/usr/bin/env Rscript

# Generate F4 ratio estimates of α (nomenclature from Patterson et al. 2012 Genetics) per chromosome separately for recently and historically admixed Amboseli animals (i.e., treat them as two separate X populations). α = the contribution of population B’ (unknown ancestral population represented by modern population B) to population X. Here, α will be an estimate of anubis admixture proportions so we'll define population B as SNPRC anubis founders.

library(admixr)

# Set data prefix so R knows where to find our eigenstrat formatted files and what they will be called
data_prefix <- "./eig.n49.CHROM"
snps <- eigenstrat(data_prefix)

# Define our admixed populations (population X) 
pops1 <- c("Amboseli_historic", "Amboseli_recent")

# Generate estimates of α for the following phylogenetic configurations (A, B, C, O): (i) hamadryas, SNPRC anubis, Mikumi yellow, macaque; (ii) Guinea, SNPRC anubis, Mikumi yellow, macaque; (iii) hamadryas, SNPRC anubis, Mikumi yellow, gelada; and (iv) Guinea, SNPRC anubis, Mikumi yellow, gelada.

# Macaque as population O/outgroup
result1 <- f4ratio(X = pops1, A = "Hamadryas", B = "SW_olive", C = "Mikumi", O = "macaque", data = snps) # phylogenetic configuation (i)
result2 <- f4ratio(X = pops1, A = "Guinea", B = "SNPRC_anubis_founder", C = "Mikumi", O = "macaque", data = snps) # phylogenetic configuation (ii)

# Gelada as population O/outgroup
result3 <- f4ratio(X = pops1, A = "Hamadryas", B = "SNPRC_anubis_founder", C = "Mikumi", O = "Gelada", data = snps) # phylogenetic configuation (iii)
result4 <- f4ratio(X = pops1, A = "Guinea", B = "SNPRC_anubis_founder", C = "Mikumi", O = "Gelada", data = snps) # phylogenetic configuation (iv)

# Combine results from all phylogenetic configurations
result_all <- rbind(result1, result2, result3, result4)
result_all$chrom <- c("CHROM") # add chromosome identifier to data frame

save(result_all, file="f4_ratio_CHROM.Rd") # save results as an R data file

print("done")
q(save="no")
