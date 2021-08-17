**Ancestry HMM (Corbett-Detig & Nielsen 2017 _PLOS Genetics_)** uses a hidden markov model and read counts per allele (read depth, extracted from a vcf file) to identify local ancestry. Our simulations reveal that it performs best on high coverage data and worse than LCLAE on low coverage individuals like most of the dataset we sample. Nevertheless, ancestry calls using Ancestry HMM are able to recapitulate all major analyses including signatures of selection against anubis introgression. 

Here we provide code to call local ancestry using Ancestry HMM which was applied to both simulated and Amboseli individuals. In the first code chunk we prepare a data file for Ancestry HMM starting with a vcf file of reference and test individuals. In part 2, we then run Ancestry HMM and create ancestry tracts from the posterior calls. The resulting files have five columns (chromosome, start, end, ancestry state, length) and can be used in place of LCLAE ancestry calls for later analyses (e.g. estimating mean anubis ancestry). 

```console 
## start with merged vcf file containing test and reference individuals
ls ./test_and_refpanel.vcf.gz 

### Prep data from vcf file, splitting each individual into two columns of read depths 
module load bcftools; bcftools query -f '%CHROM \t %POS \t %REF \t %ALT [\t %AD]\n' ./test_and_refpanel.vcf.gz  > temp
sed -i 's/\./0,0/g' temp
sed -i 's/,/\t/g' temp

module load bcftools ; bcftools query -l ./test_and_refpanel.vcf.gz  > samples.list

## get the snps to be considered
cut -f 1-2 temp > temp2.sites 

## samples start in column 5, with each individual being represented by 2 columns
## get test data (e.g. samples 1-10 in the vcf which will be columns 5-24)
cut -f 5-24 temp > temp2.testdata
## get reference data (e.g. samples 11-30 in the vcf which will be columns 25-64
cut -f 25-64 temp > temp2.refdata
head -30 samples.list | tail -20 > refsamples.list

## in R, sum allele counts for the reference panel populations 
## alternately, the mean allele frequency can be used instead of allele counts if the reference panel is of sufficiently high coverage
module load R; R
library(data.table); fread("temp2.sites", header=F) -> sites; fread("temp2.refdata", header=F) -> data; colnames(sites) <- c('chrom', 'snp')
## separate refpanel alleles by population. for us, ref1 was yellow and ref2 was anubis
ref1 <- data[,c(3:4, 11:14, 43:44, 51:56)]; ref2 <- data[,!c(3:4, 11:14, 43:44, 51:56)] 
## sum reference and alternate alleles
keep <- seq(1,ncol(ref2),2); sites$ref_ref2 <- rowSums(ref2[,..keep], na.rm=T); keep <- seq(2,ncol(ref2),2); sites$alt_ref2 <- rowSums(ref2[,..keep], na.rm=T); keep <- seq(1,ncol(ref1),2); sites$ref_ref1 <- rowSums(ref1[,..keep], na.rm=T); keep <- seq(2,ncol(ref1),2); sites$alt_ref1 <- rowSums(ref1[,..keep], na.rm=T)
## calculate distance between SNPs assuming a uniform recombination rate of 1 cM/Mb (or 1 M/1e8 bp). Express in M/bp. 
## we don't use the baboon recombination rates we call to avoid biasing downstream analyses
sites$dist <- sites$snp - c(1, sites$snp[-nrow(sites)]); sites$dist <- sites$dist * 1e-8
## export intermediary without test data. This was done because we were exceeding memory trying to do it all at once. 
write.table(, "data_no_test.txt", row.names=F, col.names=T, sep="\t", quote=F)
## add in testdata
fread("temp2.testdata") -> test; cbind(sites, test) -> sites
write.table(sites, "data.panel", row.names=F, col.names=F, sep="\t", quote=F)

## Get data for each chromosome 
for chrom in `cat 00chroms`; do mkdir $chrom; cp test.samples $chrom/ ; grep '^${chrom}' data.panel > $chrom/data.panel ; done
## for chrom 1 and 2, manually subset the rows of data.panel to avoid having other chromosomes which start with the same number 

## Thin to remove sites that are not differentiated between populations (minimum allele frequency difference of 20%)
module load R; R
library(data.table); fread("data.panel") -> d
ref1 <- d$V4/rowSums(d[,3:4]); ref2 <- d$V6/rowSums(d[,5:6])
d2 <- subset(d, abs(ref2-ref1) > 0.2)
write.table(d2, "filtered.panel", row.names=F, col.names=F, sep="\t", quote=F)
```



```console
## Prepare to run Ancestry HMM
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/Programs/conda/lib/
export PATH=$PATH:~/Programs/conda/lib/
module load Anaconda; conda --version
module load libtool; libtool --version

# Install the AncestryHMM software
# wget https://github.com/russcd/Ancestry_HMM/archive/master.zip; gunzip master.zip
# cd Ancestry_HMM/src/
# make

~/Programs/Ancestry_HMM/src/ancestry_hmm -i filtered.panel -s test.samples -a 2 0.25 0.75 -p 1 100000 0.75 -p 0 -100 0.25 --fixed
# -a 2 0.25 0.75 says there are 2 ancestral populations contributing 25% and 75% of ancestry to the population
# -p for each admixture pulse has 3 parameters referring to (i) the source population, (ii) how many generations ago they enterred the population (for negative values) or the effective population size (for positive values and the ancestral population), and (iii) the proportion of modern ancestry contributed by that pulse
# --fixed notes that the timing of admixture is set

# Create ancestry tracts from posterior calls
# array_for_get_ancestry.sh with INDIV and CHROM, which calls get_ancestryHMM_tracts.R 
for chrom in `cat 00chroms`; do sed -e s/CHROM/$chrom/g get_ancestryHMM_tracts.R > get.$chrom.R; sed -e s/CHROM/$chrom/g array_for_get_ancestry.sh > array.$chrom.sh; sbatch --array=1-445 array.$chrom.sh; rm array.$chrom.sh; done
# outputs a file for CHROM_INDIV in tracts

# Convert tracts to HMM tracts
for i in `cat 00_ambo_names.list`; do cat tracts/*${i}* > HMMtracts/${i}.txt; done 
```


#### Internal scripts

The above code calls the functions `array_for_get_ancestry.sh` and `get_ancestryHMM_tracts.R` which are included below. 

```console
#!/bin/bash
#SBATCH --get-user-env

index=${SLURM_ARRAY_TASK_ID}
chr=CHROM

name=`head -$index 00_ambo_names.list | tail -1`
echo $name
echo $index

sed -e s/INDIV/$name/g get.$chr.R > get.$name.$chr.R
chmod 777 get.$name.$chr.R
module load R
./get.$name.$chr.R
rm get.$name.$chr.R
```

```console
#!/usr/bin/env Rscript
#SBATCH --get-user-env

read.delim("~/genomes/panubis1/Panubis1.0.fa.fai", header=F) -> chroms
chr_length <- chroms$V2[chroms$V1 == "CHROM"]

library(data.table); fread("CHROM/INDIV.posterior") -> data
colnames(data)[3:5] <- c("Anu", "Het", "Yel")
data$maxs <-  apply(data[,3:5],1, max)
data$state <- NA

data$state[data$Anu == data$maxs] <- 2
data$state[data$Het == data$maxs] <- 1
data$state[data$Yel == data$maxs] <- 0

data$nxt_pos <- c(data$position[-1], chr_length)
data$nxt_state <- c(data$state[-1], data$state[nrow(data)])
d2 <- data[!data$state == data$nxt_state,]

d2[,.(chrom, position, nxt_pos, state, nxt_state)] -> blocks
###Let's add the first block from position 1 to the first break pt
  first <- as.data.frame(t(matrix(c(blocks$chrom[1], 1, blocks$position[1], blocks$state[1], blocks$state[1]))))
  colnames(first) <- colnames(blocks); rbind(first, blocks) -> blocks; rm(first)

## Start of each tract is the midpoint between position and nxt_pos, and that tract is characterized by the 'nxt_state'
blocks$start <- apply(blocks[,2:3],1,mean)
blocks$start[1] <- 1
blocks$end <- c(blocks$start[-1],chr_length)
blocks$length <- blocks$end - blocks$start
toprint <- blocks[,c(1,6,7,5,8)]
write.table(toprint, "./tracts/CHROM_INDIV.tracts", row.names=F, col.names=F, sep="\t", quote=F)
```
