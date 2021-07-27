#!/usr/bin/env Rscript

library(admixr)

# set data prefix so R knows where to find our eigenstrat formatted files
data_prefix <- "./eig.n51.CHROM"
snps <- eigenstrat(data_prefix)

# define our populations when our B populaton is anubis (we want to estimate anubis admixture in Amboseli)
pops1 <- c("Amboseli_historic", "Amboseli_recent")

result1 <- f4ratio(X = pops1, A = "Guinea", B = "SW_olive", C = "Mikumi", O = "macaque", data = snps)
result2 <- f4ratio(X = pops1, A = "Hamadryas", B = "SW_olive", C = "Mikumi", O = "macaque", data = snps)

# repeat with gelada as the outgroup instead of macaque
result3 <- f4ratio(X = pops1, A = "Guinea", B = "SW_olive", C = "Mikumi", O = "Gelada", data = snps)
result4 <- f4ratio(X = pops1, A = "Hamadryas", B = "SW_olive", C = "Mikumi", O = "Gelada", data = snps)

result_all <- rbind(result1, result2, result3, result4)
result_all$chrom <- c("CHROM")

save(result_all, file="f4_ratio_CHROM.Rd")

print("done")
q(save="no")
