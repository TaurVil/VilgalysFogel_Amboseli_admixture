# Script for evaluating whether anubis allele frequency change over time is predicted by local genomic features (Supplementary Methods Section 14.2)

# Load R libraries
library(dplyr)
library(ggplot2)
library(data.table)

# Load ancestry in 100 kb windows for all Amboseli individuals
load("./100kb_ancestry_and_features.RData")

ancestry_windows <- cbind(new_windows, window_by_individual)

# Load demographic data for Amboseli individuals with confirmed ids, including the starting and endings years they were present in the population
# Other demographic data used in other scripts are also included in this data frame but which we won't use here (their sex, entry type (B = born into a study group, O = first observed in a new study group, or I = immigrating into a study group)
demog <- read.table("amboseli.demographic.info.txt", header=T)

# Get a list of individuals in the population that were present in each year
years_list <- seq(min(demog$start_year), max(demog$end_year)) # get list of all unique years in the demographic data
max(years_list) # 2021
# Remove 2021 because we do not have data for the entire year
years_list <- years_list[1:length(years_list)-1]

for (i in 1:length(years_list)) { # run for each year in the data set
  
  tmp <- demog[demog$start_year<=years_list[i] & demog$end_year>=years_list[i],] # for each year, get individuals present in the population during that year
  
  tmp2 <- tmp[c(1)] # get individual ids 
  tmp2[,2] <- years_list[i] # add focal year to the data frame
  
  if (i==1) {tmp2 -> years}
  if (i>1) {rbind(years, tmp2) -> years}
  
}

colnames(years)[2] <- c("year")

# How many individuals are present in the population every year?
tmp <- years %>% group_by(year) %>% tally()
# Visualize the number of individuals present in the population every year
ggplot(data=tmp) +geom_col(aes(year, n), color="white") + geom_text(aes(year,n, label=n), vjust=-1, size=3) + theme_classic() + theme(text=element_text(size=12), axis.text = element_text(color="black")) #+ geom_vline(aes(xintercept=1979-0.5), linetype="dashed", color="steelblue3", size=1) # we have between 1-5 individuals through 1978 and then the sample size approximately doubles in 1979 to 11 

# We'll start our analysis in 1979 where we begin to have more than just a handful of individuals so that we can get an estimate of average anubis ancestry in the population per year
tmp <- tmp[tmp$year>=1979,]
min(tmp$n) # minimum of 11 individuals resequenced per year
max(tmp$n) # maximum of 228 individuals resequenced per year
round(mean(tmp$n),3) # mean of 135.19 individuals resequenced per year

which(years_list>=1979) # 9-50
years_list <- years_list[9:50]

# For each year, get the average anubis ancestry in the population in each 100 kb window
for (i in 1:length(years_list)) {
  
  tmp <- years[years$year==years_list[i],] # get the list of individuals to include for the focal year
  which(colnames(ancestry_windows) %in% tmp$table_s1_id) -> k
  
  tmp2 <- ancestry_windows[c(1:3,k)]
  tmp2$na_count <- rowSums(is.na(tmp2[c(4:ncol(tmp2))])) # count any NAs 
  tmp2$pop_avg_anubis_admix <- (rowSums(tmp2[c(4:(ncol(tmp2)-1))], na.rm=T))/(length(k)-tmp2$na_count) # get the average anubis ancestry for each 100 kb genomic window (accounting for any individuals with NAs)
  tmp2$indiv_count <- length(k) # add the number of individuals included in the focal year in the population (the number of individuals actually used for estimating average anubis ancestry at each site would be indiv_count-na_count)
  
  tmp3 <- tmp2[names(tmp2) %in% c("chr", "start", "end", "na_count", "pop_avg_anubis_admix", "indiv_count")] # only keep the columns we need
  tmp3[,7] <- years_list[i] # add the focal year
  
  if (i==1) {tmp3 -> ancestry_year} # if this is the first year, create a new data frame called ancestry_year
  if (i>1) {rbind(ancestry_year, tmp3) -> ancestry_year} # if this is not the first year, add the focal year's data to the existing data frame ancestry_year
  
  
  print(i)
  
}

colnames(ancestry_year)[7] <- c("year")

# Get the middle position in each window
ancestry_year$pos <- ((ancestry_year$end-ancestry_year$start)/2)+ancestry_year$start
# Make a unique id per window
ancestry_year$id <- paste(ancestry_year$chr, ancestry_year$pos, sep="_")

# For each genomic window, calculate the change in anubis ancestry over time (get the beta for the year effect in a linear model of mean anubis ancestry per year ~ year) 
# Get a list of the unique genomic windows we'll be evaluating
ancestry2 <- ancestry_windows %>% select(chr, start, end)

# Same as with the ancestry_year data frame - get the middle position in each window
ancestry2$pos <- ((ancestry2$end-ancestry2$start)/2)+ancestry2$start
# Make a unique id per window
ancestry2$id <- paste(ancestry2$chr, ancestry2$pos, sep="_")

ancestry_change <- ancestry2 

for (i in 1:nrow(ancestry2)) { # for each genomic window
  
  tmp <- ancestry_year[ancestry_year$id==ancestry2$id[i],]
  ancestry_change[i,6] <- lm(pop_avg_anubis_admix~year, data=tmp)$coeff[2] # grab the beta for year
  tmp2 <- tmp[tmp$year==min(tmp$year),]
  ancestry_change[i,7] <- tmp2$pop_avg_anubis_admix # grab the average anubis ancestry in the population in the starting year of the analysis (1979)
  ancestry_change[i,8] <- nrow(tmp) # number of years included in the model (should be 42 years for every genomic window since we have no NAs for average anubis ancestry in the population in the ancestry_year data frame)
  
  if (i/500 == round(i/500)) {print(i)}
}

colnames(ancestry_change)[6:8] <- c("annual_anubis_change_beta", "starting_pop_avg_anubis_admix_1979", "year_count")

min(ancestry_change$year_count) # all genomic windows used 42 years (what we expect)

# Add genomic features to the ancestry_change data frame so we can evaluate whether the change in anubis ancestry over time is not only predicted by the starting anubis frequency in 1979, but also FST, mean recombination rate, and B values
# 1. FST - to get these estimates, use the script 4bFST_SNPRCref.sh in "VilgalysFogel_Amboseli_admixture/02_local-ancestry-calling/02.3.pedigree_trios/" but substituting 100000 for window size --fst-window-size in place of 35000.
fst <- read.table("fst_masked_unmerged_SNPRCref_100kbwin_500bpstep.windowed.weir.fst", header=T) 
# Remove mean FST as we will use weighted Fst (although results are qualitatively similar whether you use one measure or the other)
fst <- fst[c(1:5)] 
# Get the center position for each 100 kb window
fst$pos <- ((fst$BIN_END-fst$BIN_START+1)/2)+fst$BIN_START-1
# Give each site a unique id
fst$id <- paste(fst$CHROM, fst$pos, sep="_")
# Remove chromosome, bin start, bin end, and position columns as we'll only need the id column for merging with the ancestry_change data frame
fst <- fst[c(4,5,7)]

# Merge FST estimates with ancestry_change data frame
all_data <- merge(ancestry_change, fst, by=c("id"), all.x=T)
summary(all_data) # we have three genomic windows with no FST estimates
subset(all_data, is.na(N_VARIANTS)) # looks like it is because there are no variants - but how do we have ancestry in those three windows then? if you check the ancestry tracts file ("amboseli_LCLAE_tracts.txt"), many individuals have a single tract that spans the entire window. Otherwise, looks like they have two tracts flanking the start and end.
# chr2_84375000 chr2 84250000 84500000 84375000 
# nrow(subset(tracts, chrom=="chr2" & start<84250000 & end>84500000)) # 323 individuals have a single tract that spans the entire window
# chr2_84625000 chr2 84500000 84750000 84625000
# nrow(subset(tracts, chrom=="chr2" & start<84500000 & end>84750000)) # 338
# chr4_34625000 chr4 34500000 34750000 34625000 
# nrow(subset(tracts, chrom=="chr4" & start<34500000 & end>34750000)) # 245

# 2. Recombination rate and B values - in the new_windows data frame but need to remove windows with recombination rates > 100 x median rcr (denoted by "max_rcr" value in RData file)
new_windows2 <- new_windows[new_windows$rcr<max_rcr,]

# Merge recombination rate and B values with all_data data frame
all_data <- merge(all_data, new_windows2, by=c("chr", "start", "end"))

# Run linear model evaluating whether the change in anubis ancestry over time is predicted by the starting anubis frequency in 1979, FST, log10 transformed mean recombination rate, and B values
summary(lm(annual_anubis_change_beta~starting_pop_avg_anubis_admix_1979+log10(rcr)+WEIGHTED_FST+B, data=all_data))

# Write model output to csv file
write.csv(summary(lm(annual_anubis_change_beta~starting_pop_avg_anubis_admix_1979+log10(rcr)+WEIGHTED_FST+B, data=all_data))$coef, file="tableS6_results.csv")

