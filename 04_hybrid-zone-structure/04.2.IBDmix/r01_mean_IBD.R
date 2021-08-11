library(data.table); library(ggplot2)
min_length <- 50000; min_lod <- 10

read.delim("~/genomes/panubis1/chromInfo.txt", header=F) -> chroms; colnames(chroms)[1:2] <- c("name", "length") # Read in chromosome file 

#### read in file of test and source individuals ###############
read.delim("./IBDmix/00_yellow_sources.list", header=F) -> yel_source; read.delim("./IBDmix/00_yel.list", header=F) -> yel
read.delim("./IBDmix/00_anubis_sources.list", header=F) -> anu_source; read.delim("./IBDmix/00_anu.list", header=F) -> anu

# For each source individual, and each chromosome within a source individual, read in the resulting data file. These files are available upon request or can be generated from the above scripts. 
IBD_yellow <- NULL; for (i in 1:nrow(yel_source)) { tmp <- NULL; for (chrom in 1:20) {
    name=paste("~/Baboon/Paper1b_demographicinference/IBDmix/IBDmix_results/yellow.relative_to_",yel_source[i,1],".",chrom,".txt",sep="")
    fread(name) -> data; data$length <- data$end - data$start; data$source <- yel_source[i,1]
    data <- subset(data, data$length >= min_length & data$slod >= min_lod); rbind(tmp,data) -> tmp }
  rbind(IBD_yellow,tmp) -> IBD_yellow; rm(tmp); print(i) }; rm(i, chrom, data, name)
IBD_yellow <- IBD_yellow[order(IBD_yellow$ID,IBD_yellow$source,IBD_yellow$chrom,IBD_yellow$start),]

IBD_anubis <- NULL; for (i in 1:nrow(anu_source)) { tmp <- NULL;  for (chrom in 1:20) {
    name=paste("~/Baboon/Paper1b_demographicinference/IBDmix/IBDmix_results/anubis.relative_to_",anu_source[i,1],".",chrom,".txt",sep="")
    fread(name) -> data; data$length <- data$end - data$start; data$source <- anu_source[i,1]
    data <- subset(data, data$length >= min_length & data$slod >= min_lod); rbind(tmp,data) -> tmp; rm(data) }
  rbind(IBD_anubis,tmp) -> IBD_anubis; rm(tmp); print(i) }; rm(i, chrom, name)
IBD_anubis <- IBD_anubis[order(IBD_anubis$ID,IBD_anubis$source,IBD_anubis$chrom,IBD_anubis$start),]

# For each individual, calculate the proportion ancestry for each source individual
colnames(yel) <- colnames(yel_source) <- "name"
for (i in yel$name) { for (j in yel_source$name) {
    subset(IBD_yellow, IBD_yellow$ID == i & IBD_yellow$source == j) -> tmp
    yel[[paste(j)]][yel$name == i] <- sum(tmp$length)/sum(chroms$length) } }; rm(tmp, i,j)
colnames(anu) <- colnames(anu_source) <- "name"
for (i in anu$name) { for (j in anu_source$name) {
    subset(IBD_anubis, IBD_anubis$ID == i & IBD_anubis$source == j) -> tmp
    anu[[paste(j)]][anu$name == i] <- sum(tmp$length)/sum(chroms$length) } }; rm(tmp, i,j)

# Get mean IBD between each yellow baboon and possible source populations
yel2 <- yel[,1:2]; colnames(yel2) <- c("individual", "population")
yel2$population[c(1:10,18)] <- 'Mikumi'; yel2$population[11:17] <- 'SNPRCyellow'
yel2$SNPRCanubis <- rowMeans(yel[,c(3:25,28:29)])
yel2$guinea <- rowMeans(yel[,35:36])
yel2$hamadryas <- rowMeans(yel[,30:31])
yel2$aberdares <- rowMeans(yel[,26:27])

# Report numbers for mean shared ancestry between yellow and anubis baboons, split by the yellow population. We use SNPRCanubis rather than anubis from the Aberdares because at least one indivdual in the Aberdares is predicted to contain considerable yellow ancestry (Rogers et al. 2019). 
mean(unlist(yel[yel2$population == "Mikumi", colnames(yel2) == "SNPRCanubis"]), na.rm=T)
sd(unlist(yel[yel2$population == "Mikumi", colnames(yel2) == "SNPRCanubis"]), na.rm=T)
mean((yel2[yel2$population == "SNPRCyellow", colnames(yel2) == "SNPRCanubis"]), na.rm=T)
sd((yel2[yel2$population == "SNPRCyellow", colnames(yel2) == "SNPRCanubis"]), na.rm=T)

# Get mean IBD between each anubis baboon and possible source populations
anu2 <- anu[,1:2]; colnames(anu2) <- c("individual", "population")
anu2$population <- "SNPRCanubis"; anu2$population[25:26] <- 'Aberdares'
anu2$SNPRCyellow <- rowMeans(anu[,12:18])
anu2$mikumi <- rowMeans(anu[,c(2:11,19)])
anu2$kinda <- rowMeans(anu[,c(22:24)])
anu2$chacma <- rowMeans(anu[,27:28])

# Report numbers for anubis individual estimated to contain recent anubis ancestry by Rogers et al. 2019, and for all other anubis baboons. 
mean(anu2[anu2$individual == "panu_30877", colnames(anu2) == "mikumi"])
mean(anu2[!anu2$individual == "panu_30877", colnames(anu2) == "mikumi"])
sd(anu2[!anu2$individual == "panu_30877", colnames(anu2) == "mikumi"])

write.table(yel2, "./ibdmix_yellow_estimates.txt", row.names=T, col.names=T, quote=F, sep="\t") # used for figure 1C and in the SI tables
write.table(anu2, "./ibdmix_anubis_estimates.txt", row.names=T, col.names=T, quote=F, sep="\t") # used in the SI tables
