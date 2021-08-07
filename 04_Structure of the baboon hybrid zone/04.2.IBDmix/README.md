## Estimate identity by descent (IBD) between anubis and yellow baboons. 

Results from LCLAE suggest some populations of yellow baboons previously thought to be unadmixed, specifically baboons used to found the SNPRC colony, may actual contain considerably anubis ancestry. However, these analyses rely upon reference panels of supposedly unadmixed individuals and admixture may extend further into the yellow baboon reference panel. Therefore, we used IBDmix (Chen et al. 2020, https://doi.org/10.1016/j.cell.2020.01.012) to estimate the amount of shared ancestry between yellow and anubis baboon individuals of different populations. IBDmix is a reference-free approach originally developed to detect archaic Neanderthal ancestry in modern African populations, meaning it does not assume any individuals are unadmixed. Specifically, IBDmix identifies sections of the genome which are IBD between individuals within a test population and a single individual representing the source of introgression. Because IBD is bidirectional (meaning it is unable to differentiate between anubis ancestry in yellow baboons and yellow ancestry in anubis baboons), we test all possible pairs of yellow and anubis baboons. We look for signatures of IBD that are shared across many potential source indivdiuals, reasoning that introgressed ancestry is unlikely to rise to high allele frequencies outside of the commonly recognized hybrid zone. 

A second limitation to IBD-based inference of introgression is that a baseline level of IBD is expected due to incomplete lineage sorting. To understand the expected rates of IBD in yellow an anubis baboons, we also apply IBDmix to estimate shared ancestry between yellow baboons with Guinea and hamadryas baboons, and between anubis baboons with yellow and Kinda baboons. These six species diverged approximately 1.5 million years ago into two clades: the northern clade of Guinea, anubis, and hamadryas baboons, and the southern clade of yellow, Kinda, and chacma baboons. As these other species are unlikely to have hybridized outside of their clade (due to strong geographic separation), estimates of IBD between different species serve as a proxy for the expected amount of IBD under incomplete lineage sorting alone without admixture. 


### Installation of IBDmix
```console 
# Local installation of up to date cmake 
./configure  --prefix=~/Programs/cmake/
gmake; gmake install 
# add cmake to the path
export PATH=$PATH:~/Programs/cmake/bin/

# install ibdmix: https://github.com/PrincetonUniversity/IBDmix
git clone https://github.com/PrincetonUniversity/IBDmix.git
cd IBDmix; mkdir build; cd build; cmake ..; cmake --build

# generates the `generate_gt` and `ibdmix` executables
~/Programs/IBDmix/IBDmix-master/build/src/
```

### Run IBDmix
We run IBDmix in parallel across chromosomes. For each run, we first grab genotype calls from previously called files (see Section 1) for yellow, anubis, and other source individuals. Genotype calls are then integrated for each data set using `generate_gt`, a function within IBDmix. IBDmix is then run sequentially for each source file and the results are exported as a series of text files.  

```console 
## have "source" files which contain the source individuals we want to iterate calling local ancestry over: 00_yellow_sources.list, 00_anubis_sources.list
## for each chromosome, run "run_IBDmix.sh"
sbatch --array=1-20 --mem=16G run_IBDmix.sh
```

### Estimate the mean proportion of the genome IBD between each test individual and source population
Following suggestions in Chen et al., we filter for tracts of IBD at least 50kb in length and with a LOD score greater than 10 in order to estimate the proportion of the genome shared between yellow and anubis baboons. 



```console
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

#### For each individual, calculate the proportion ancestry for each source individual ############
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
# yel2$sd_sw <- apply(yel[,c(3:25,28:29)], 1, sd)/sqrt(length(c(3:25,28:29)))
# yel2$sd_aberdare <- apply(yel[,26:27], 1, sd)/sqrt(length(26:27))
# yel2$sd_ham <- apply(yel[,30:31], 1, sd)/sqrt(length(30:31))
# yel2$sd_gui <- apply(yel[,35:36], 1, sd)/sqrt(length(35:36))

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
```

### I

Next we sought to identify the overlap between shared ancestry identified using LCLAE and IBDmix. 

.  When looking for overlap between multiple source individuals we applied more liberal thresholds for calling a region of the genome IBD, specifically tracts longer than 1000bp and LOD score greater than 4. While these criteria likely retain false positives, 


### For comparison to the SNPRC "yellow" founders, estimate IBD between Amboseli and all baboon species
We repeated the above process for nine high coverage Amboseli baboons to serve as a positive control of what , using all baboon species profiled as possible sources of ancestry. The files to run IBDmix are included here as `00_amboseli.list`, `00_amboseli_sources.list`, and `run_IBDmix_amboseli.sh`. Allele sharing between Amboseli baboons and each population was estimating using the follow R code. Unlike above, we don't collect where in the genome these tracts are located but rather simply estimate a mean proportion of the genome which is IBD between each Amboseli baboon and other baboon species. 

```console
library(data.table); library(ggplot2)
min_length <- 50000; min_lod <- 10 # thresholds for minimum length and LOD score

# Read in chromosome file, with the chromosome names and length
read.delim("~/genomes/panubis1/chromInfo.txt", header=F) -> chroms; colnames(chroms)[1:2] <- c("name", "length")

# Read in lists of test (Amboseli) and source (all other baboon) individuals
read.delim("./IBDmix/00_amboseli_sources.txt", header=F) -> sources
read.delim("./IBDmix/00_amboseli.list", header=F) -> samples

# For each source individual, and each chromosome within a source individual, read in the resulting data file. These files are available upon request or can be generated from the above scripts. 
IBD_ambo <- NULL; for (i in 1:nrow(sources)) {
  tmp <- NULL
  for (chrom in 1:20) {
    name=paste("./IBDmix_amboseli/amboseli.relative_to_",sources[i,1],".",chrom,".txt",sep="")
    fread(name) -> data; data$length <- data$end - data$start; data$source <- sources[i,1]
    data <- subset(data, data$length >= min_length & data$slod >= min_lod)
    rbind(tmp,data) -> tmp
  }
  rbind(IBD_ambo,tmp) -> IBD_ambo; rm(tmp); print(i)
}; rm(i, chrom, data, name)
IBD_ambo <- IBD_ambo[order(IBD_ambo$ID,IBD_ambo$source,IBD_ambo$chrom,IBD_ambo$start),]

# For each individual, calculate the proportion ancestry for each source individual
colnames(samples) <- colnames(sources) <- "name"
for (i in samples$name) {
  for (j in sources$name) {
    subset(IBD_ambo, IBD_ambo$ID == i & IBD_ambo$source == j) -> tmp
    samples[[paste(j)]][samples$name == i] <- sum(tmp$length)/sum(chroms$length)
  }
}; rm(tmp, i,j)

# Summarize output by Amboseli individual and source population
amboseli <- samples[,1:2]
colnames(amboseli) <- c("individual", "population")
amboseli$population <- "Amboseli"

amboseli$mikumi <- rowMeans(samples[,c(30:39, 47)])
amboseli$hamadryas <- rowMeans(samples[,c(48:49)])
amboseli$kinda <- rowMeans(samples[,c(50:52)])
amboseli$guinea <- rowMeans(samples[,c(53:54)])
amboseli$chacma <- rowMeans(samples[,c(55:56)])
amboseli$aberdare <- rowMeans(samples[,c(26:27)])
amboseli$SNPRCyellow <- rowMeans(samples[,c(40:46)])
amboseli$SNPRCanubis <- rowMeans(samples[,c(2:25,28:29)])

write.table(amboseli, "./ibdmix_amboseli_estimates.txt", row.names=F, col.names=T, sep="\t", quote=F)

```
