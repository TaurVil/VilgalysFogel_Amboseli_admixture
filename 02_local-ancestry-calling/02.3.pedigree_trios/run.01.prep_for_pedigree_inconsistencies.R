# First script for evaluating the quality of our local ancestry calls from LCLAE using known parent-offspring trios from the Amboseli population (see Supplementary Methods 7.3)

# Load R libraries needed
library(dplyr)

# Load text file containing pedigree trios and their associated genomic data (genome-wide anubis ancestry and coverage)
trios <- read.table("pedigree_trio_info.txt", header=T)

# Load ancestry tracts for Amboseli individuals
tracts <- read.table("amboseli_LCLAE_tracts.txt", header=T)

# Remove tracts <1 kb in length
tracts <- tracts[tracts$length2>1000,]

# Retain tracts only for individuals in the pedigree trios
# First, get list of unique individuals across all pedigree trios
tmp <- trios[c(1)] # kid's id
tmp2 <- trios[c(2)] # mom's id
tmp3 <- trios[c(3)] # dad's id
# Rename column names so we can combine all three into one data frame
colnames(tmp) <- c("indiv")
colnames(tmp2) <- c("indiv")
colnames(tmp3) <- c("indiv")
tmp4 <- rbind(tmp, tmp2, tmp3)
# Get list of unique individuals across all moms, dads, and kids
indiv_list <- tmp4 %>% distinct(indiv)
nrow(indiv_list) # 181 unique individuals

# Merge tracts with unique individuals across the pedigree trios so we only retain tracts from individuals we'll work with in this analysis
tmp2 <- merge(indiv_list, tracts, by.x=c("indiv"), by.y = c("table_s1_id"))

# Check to make we have tract data for all 181 individuals
nrow(distinct(tmp2, indiv))==nrow(indiv_list) # TRUE

# Grab only columns we'll need (indiv, chrom, start of the tract, end of the tract, ancestry state of the tract)
tracts2 <- tmp2[c(1:5)]

# We also need to get the positions along the genome we'd like to evaluate for pedigree consistencies/inconsistencies
# Load chromosomes sizes from the Panubis1 genome
chrom_sizes <- read.table("chromsizes_Panubis1_final.txt", header=FALSE)
colnames(chrom_sizes) <- c("chrom", "length")

# Remove chrX and chrY
chrom_sizes <- chrom_sizes[chrom_sizes$chrom!="chrY" & chrom_sizes$chrom!="chrX",]

# Exclude the first and last 50 kb of each chromosome (and account for the fact that we'll use a window size of 35 kb)
# Define window size
wind=35000 # use 35 kb as our window size because this is the same window used for calling ancestry in LCLAE
chrom_sizes$start_chrom <- 50000+(wind/2)
chrom_sizes$end_chrom <- chrom_sizes$length-(50000+(wind/2))

# Grab only the columns we'll need (i.e., only the chromosome name, start of chromosome, and end of chromosome so exclude the chromosome length column)
chrom_sizes <- chrom_sizes[c(1,3:4)]

# Get a list of the chromosomes we'll look at (chromosomes 1-20)
chrom_list <- chrom_sizes[c(1)]

# Generate list of positions along the genome where positions are 35 kb apart along each chromosome starting 50 kb from the beginning of the chromosome and ending at least 50 kb from the end of the chromosome)
for (i in 1:nrow(chrom_list)) {
  
  tmp <- chrom_sizes[i,]
  tmp2 <- as.data.frame(seq(from=tmp$start_chrom, to=tmp$end_chrom, by=wind))
  tmp2[,2] <- tmp[,1]
  tmp3 <- tmp2[c(2,1)]
  colnames(tmp3) <- c("chrom", "pos")
  
  if (i==1) {tmp3 -> positions}
  if (i>1) {rbind(positions, tmp3) -> positions}
  print(i)
  
}

# How many positions will we look at across the genome?
nrow(positions) # 73975 

# Check to confirm that 73975 is the correct number of total positions 
sum(floor(((chrom_sizes$end_chrom-chrom_sizes$start_chrom)/wind)+1))==nrow(positions) # TRUE
# Now we have all of the positions where we want evaluate consistences of ancestry calls in the pedigree trios

# Save dataframes for upload to the computing cluster for parallelization 
save(indiv_list, tracts2, positions, trios, file="local_ancestry_pedigree_trios_maskedSNPRCref.Rd")
