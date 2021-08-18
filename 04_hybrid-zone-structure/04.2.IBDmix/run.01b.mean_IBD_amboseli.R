library(data.table); library(ggplot2)
min_length <- 50000; min_lod <- 10 # thresholds for minimum length and LOD score

# Read in chromosome file, with the chromosome names and length. This is the first two columns of the Panubis1 genome idex
read.delim("~/genomes/panubis1/chromInfo.txt", header=F) -> chroms; colnames(chroms)[1:2] <- c("name", "length")

# Read in lists of test (Amboseli) and source (all other baboon) individuals
read.delim("./DATA/00_amboseli_sources.list", header=F) -> sources; read.delim("./DATA/00_amboseli.list", header=F) -> samples

# For each source individual, and each chromosome within a source individual, read in the resulting data file. These files are available upon request or can be generated from the above scripts. 
IBD_ambo <- NULL; for (i in 1:nrow(sources)) { tmp <- NULL; for (chrom in 1:20) {
    name=paste("./IBDmix_by_chrom/amboseli.relative_to_",sources[i,1],".",chrom,".txt",sep="")
    fread(name) -> data; data$length <- data$end - data$start; data$source <- sources[i,1]
    data <- subset(data, data$length >= min_length & data$slod >= min_lod); rbind(tmp,data) -> tmp }
  rbind(IBD_ambo,tmp) -> IBD_ambo; rm(tmp); print(i) }; rm(i, chrom, data, name)
IBD_ambo <- IBD_ambo[order(IBD_ambo$ID,IBD_ambo$source,IBD_ambo$chrom,IBD_ambo$start),]

# For each individual, calculate percent of the genome IBD with each source individual
colnames(samples) <- colnames(sources) <- "name"
for (i in samples$name) { for (j in sources$name) { subset(IBD_ambo, IBD_ambo$ID == i & IBD_ambo$source == j) -> tmp
    samples[[paste(j)]][samples$name == i] <- sum(tmp$length)/sum(chroms$length) } }; rm(tmp, i,j)

# Summarize output by Amboseli individual and source population
amboseli <- samples[,1:2]; colnames(amboseli) <- c("individual", "population"); amboseli$population <- "Amboseli"

amboseli$mikumi <- rowMeans(samples[,c(30:39, 47)])
amboseli$hamadryas <- rowMeans(samples[,c(48:49)])
amboseli$kinda <- rowMeans(samples[,c(50:52)])
amboseli$guinea <- rowMeans(samples[,c(53:54)])
amboseli$chacma <- rowMeans(samples[,c(55:56)])
amboseli$aberdares <- rowMeans(samples[,c(26:27)])
amboseli$SNPRCyellow <- rowMeans(samples[,c(40:46)])
amboseli$SNPRCanubis <- rowMeans(samples[,c(2:25,28:29)])

write.table(amboseli, "./RESULTS/ibdmix_amboseli_estimates.txt", row.names=F, col.names=T, sep="\t", quote=F)
