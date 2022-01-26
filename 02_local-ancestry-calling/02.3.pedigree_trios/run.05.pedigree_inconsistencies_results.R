# Final script for evaluating the quality of our local ancestry calls from LCLAE using known parent-offspring trios from the Amboseli population (see Supplementary Methods 7.3)

# Load R libraries
library(ggplot2)
library(data.table)
library(lmerTest)

# Load pedigree trio data, ancestry state inconsistency results, and covariate data for the SNPRC reference panel
load("local_ancestry_pedigree_trios_maskedSNPRCref.Rd", verbose=T) # Pedigree trio data
#Loading objects:
# indiv_list
# tracts2
# positions
# trios
load("inconsistent_local_ancestry_calls_pedigree_trios_maskedSNPRCref.Rd", verbose=T) # Ancestry state inconsistency results 

#Loading objects:
#  inconsistent - 1 is inconsistent, 0 is not inconsistent and trio is usable, NA is trio is not usable due to not having a call for one of the individuals (boundary, not in tract)
#  inconsistent_error - if the pedigree trio is inconsistent and it is a case where there are two ways of trio being inconsistent (i.e., when both parents are homozygous), get what erroneous ancestry is being called
#  call_records - ancestry state calls for all members of the pedigree trio
aims <- read.table("all.for_pedigree_trios_AIM_count_maskedSNPRCref_pedigree_trios.txt", header=F) # Ancestry informative marker (AIM) counts
fst <- read.table("fst_masked_unmerged_SWref_35kbwin_500bpstep.windowed.weir.fst", header=T) # FST
load("for_pedigree_trios_250kbwin_recombination.Rd", verbose=T) # Mean recombination rate
#Loading objects:
# rcr

#############################################################################################################################
# Ancestry state inconsistencies/consistencies across pedigree trios.
#############################################################################################################################

# Give each site a unique id
inconsistent$id <- paste(inconsistent$chrom, inconsistent$pos, sep="_")

# Get the proportion of inconsistent ancestry calls across pedigree trios per site
inconsistent$total_inconsistent <- rowSums(inconsistent[c(3:(nrow(trios)+2))], na.rm=T)
inconsistent$total_NA <- rowSums(is.na(inconsistent[c(3:(nrow(trios)+2))]))
inconsistent$total_consistent <- nrow(trios)-(inconsistent$total_inconsistent+inconsistent$total_NA)
inconsistent$prop_inconsistent <- inconsistent$total_inconsistent/(inconsistent$total_inconsistent+inconsistent$total_consistent)
inconsistent$prop_NA <- inconsistent$total_NA/nrow(trios)
inconsistent$prop_consistent <- 1-inconsistent$prop_inconsistent # proportion of consistent ancestry calls (another way of thinking about these results)

# Check to make sure everything was calculated correctly
# First, sum total_inconsistent, total_consistent, and total_NA which should equal the number of pedigree trios total (=92) for all rows in the data set i.e., our total number of positions (=73,975)
sum((inconsistent$total_inconsistent+inconsistent$total_consistent+inconsistent$total_NA)==nrow(trios))==nrow(positions) # TRUE

# Second, check an individual site
tmp <- inconsistent[1,(3:(nrow(trios)+2))] # check the site in the first row - grab inconsistent data for all pedigree trios at that site
tmp2 <- as.data.frame(t(tmp))
colnames(tmp2) <- c("calls")

nrow(subset(tmp2, calls==1))==inconsistent$total_inconsistent[1] # TRUE
nrow(subset(tmp2, is.na(calls)))==inconsistent$total_NA[1] # TRUE
nrow(subset(tmp2, calls==0))==inconsistent$total_consistent[1] # TRUE

round(median(inconsistent$prop_inconsistent),3) # 0.076
round(sd(inconsistent$prop_inconsistent),3) # 0.062

round(median(inconsistent$prop_consistent),3) # 0.924
round(sd(inconsistent$prop_consistent),3) # 0.062

#############################################################################################################################
# Covariates to include in model predicting a site's propotion of ancestry state inconsistencies/consistencies across pedigree trios.
#############################################################################################################################

### 1. Ancestry informative markers
# Let's look at the number of ancestry informative markers count in 35 kb windows surrounding each focal position
colnames(aims) <- c("chrom", "pos", "start_window", "end_window", "AIM_count", "total_calls_population", "mean_calls_per_AIM", "sd_calls")

# Give each site a unique id (for merging with the other data sets)
aims$id <- paste(aims$chrom, aims$pos, sep="_")

# Remove start and end of window columns which we will not need
aims <- aims[c(1:2,5:9)]

# Merge ancestry informative marker data with pedigree inconsistency data frame
inconsistent2 <- merge(inconsistent, aims, by=c("chrom", "pos", "id"))

# Generally, how many unique ancestry informative markers do we have in 35 kb windows?
round(summary(inconsistent2$AIM_count), 2)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.0   109.0   149.0   155.6   194.0  1514.0 

# Plot a histogram of the number of unique ancestry informative markers across 35 kb windows
ggplot(inconsistent2) + geom_histogram(aes(AIM_count), color="white", alpha=0.9, bins=200) + scale_x_continuous(name="# of unique ancestry informative markers\nin 35 kb windows across the genome") + theme_classic() + theme(text=element_text(size=20)) 

### 2. FST
# Remove mean FST as we will use weighted Fst
fst <- fst[c(1:5)]
# Get the center position for each 35 kb window
fst$pos <- ((fst$BIN_END-fst$BIN_START+1)/2)+fst$BIN_START-1
# Remove start and end of window columns which we will not need
fst <- fst[c(1,4:6)]
# Give each site a unique id (for merging with the other datas sets)
fst$id <- paste(fst$CHROM, fst$pos, sep="_")

# Merge FST data with other results (this will only retain the FST windows in our analysis)
inconsistent3 <- merge(inconsistent2, fst, by.x=c("chrom", "pos", "id"), by.y=c("CHROM", "pos", "id"), all.x=T)

summary(inconsistent3$WEIGHTED_FST) # we have 600 windows with no FST values (NA's) - why?
subset(inconsistent3, is.na(WEIGHTED_FST))$N_VARIANTS # these windows do not have any variants (all NA's) which makes sense why we don't have FST values for them
sum(subset(inconsistent3, is.na(WEIGHTED_FST))$AIM_count) # 0 - these windows do not have any AIMs (all 0s) which we would expect if there are no variants

### 3. Recombination Rate
# Note that the mean recombination rate was calculated using 250 kb windows (even though our focal sites are 35 kb apart - a 35 kb window size for mean recombination rate would generate too noisy of an estimate)
# Remove any windows with 100x the median recombination rate across all 250 kb genomic windows
max_rcr <- 100*median(rcr$n24_anubis) #0.06977849
(nrow(rcr[rcr$n24_anubis > max_rcr,])/nrow(rcr))*100 # 0.07% (n=53) of our windows will have to be removed  

# Also remove windows where we do not get recombinations rates for the entire 250 kb genomic window
(nrow(subset(rcr, length_used<250000))/nrow(rcr))*100 # 0.10% (n=77) of our windows will have to be removed 

rcr <- rcr[rcr$n24_anubis<max_rcr,]
rcr <- rcr[rcr$length_used==250000,]
nrow(rcr) # left with 73,847 after applying recombination rate filters

# new maximum recombination rate
max(rcr$n24_anubis) #0.06969489
# length of the window used for calculating recombination rates
summary(rcr$length_used) # all 250000

# Give each site a unique id (for merging with the other datas sets)
rcr$id <- paste(rcr$chrom, rcr$pos, sep="_")

# Remove start and end of window columns as well as the length_used column which we will not need
rcr <- rcr[c(1:2,5,7)]

inconsistent_final <- merge(inconsistent3, rcr, by=c("chrom", "pos", "id"), all.x=T)

inconsistent_maskedSNPRCref <- inconsistent_final

rm(fst, aims, inconsistent, inconsistent_error, inconsistent2, inconsistent3, calls_record)

# How do the pedigree inconsistency results compare when using the Wall et al. reference panel?
# Mostly repeat lines 30-131 for the Wall et al. data (can skip including AIMs, FST, and rcr which are commented out - if you want to include these covariates, can use the same rcr calls so we can just merge the Wall et al. data with the rcr data frames)
load("inconsistent_local_ancestry_calls_pedigree_trios_unmaskedWallref.Rd", verbose=T) # Ancestry state inconsistency results 
#fst <- read.table("fst_unmasked_unmerged_Wallref_35kbwin_500bpstep.windowed.weir.fst", header=T)
#aims <- read.table("all.for_pedigree_trios_AIM_count_unmaskedWallref_pedigree_trios.txt", header=F) # Ancestry informative marker (AIM) counts

# Give each site a unique id
inconsistent$id <- paste(inconsistent$chrom, inconsistent$pos, sep="_")

# Get the proportion of inconsistent ancestry calls across pedigree trios per site
inconsistent$total_inconsistent <- rowSums(inconsistent[c(3:(nrow(trios)+2))], na.rm=T)
inconsistent$total_NA <- rowSums(is.na(inconsistent[c(3:(nrow(trios)+2))]))
inconsistent$total_consistent <- nrow(trios)-(inconsistent$total_inconsistent+inconsistent$total_NA)
inconsistent$prop_inconsistent <- inconsistent$total_inconsistent/(inconsistent$total_inconsistent+inconsistent$total_consistent)
inconsistent$prop_NA <- inconsistent$total_NA/nrow(trios)
inconsistent$prop_consistent <- 1-inconsistent$prop_inconsistent # proportion of consistent ancestry calls (another way of thinking about these results)

# Check to make sure everything was calculated correctly
# First, sum total_inconsistent, total_consistent, and total_NA which should equal the number of pedigree trios total (=114) for all rows in the data set i.e., our total number of positions (=73,975)
sum((inconsistent$total_inconsistent+inconsistent$total_consistent+inconsistent$total_NA)==nrow(trios))==nrow(positions) # TRUE

# Second, check an individual site
tmp <- inconsistent[1,(3:(nrow(trios)+2))] # check the site in the first row - grab inconsistent data for all pedigree trios at that site
tmp2 <- as.data.frame(t(tmp))
colnames(tmp2) <- c("calls")

nrow(subset(tmp2, calls==1))==inconsistent$total_inconsistent[1] # TRUE
nrow(subset(tmp2, is.na(calls)))==inconsistent$total_NA[1] # TRUE
nrow(subset(tmp2, calls==0))==inconsistent$total_consistent[1] # TRUE

round(median(inconsistent$prop_inconsistent),3) # 0.283
round(sd(inconsistent$prop_inconsistent),3) # 0.099

round(median(inconsistent$prop_consistent),3) # 0.717
round(sd(inconsistent$prop_consistent),3) # 0.099

### 1. Ancestry informative markers
# Let's look at the number of ancestry informative markers count in 35 kb windows surrounding each focal position
#colnames(aims) <- c("chrom", "pos", "start_window", "end_window", "AIM_count", "total_calls_population", "mean_calls_per_AIM", "sd_calls")

# Give each site a unique id (for merging with the other data sets)
#aims$id <- paste(aims$chrom, aims$pos, sep="_")

# Remove start and end of window columns which we will not need
#aims <- aims[c(1:2,5:9)]

# Merge ancestry informative marker data with pedigree inconsistency data frame
#inconsistent2 <- merge(inconsistent, aims, by=c("chrom", "pos", "id"))

# Generally, how many unique ancestry informative markers do we have in 35 kb windows?
#round(summary(inconsistent2$AIM_count), 2)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.0    80.0   110.0   116.1   146.0  1354.0 

# Plot a histogram of the number of unique ancestry informative markers across 35 kb windows
#ggplot(inconsistent2) + geom_histogram(aes(AIM_count), color="white", alpha=0.9, bins=200) + scale_x_continuous(name="# of unique ancestry informative markers\nin 35 kb windows across the genome") + theme_classic() + theme(text=element_text(size=20)) 

### 2. FST
# Remove mean FST as we will use weighted Fst
#fst <- fst[c(1:5)]
# Get the center position for each 35 kb window
#fst$pos <- ((fst$BIN_END-fst$BIN_START+1)/2)+fst$BIN_START-1
# Remove start and end of window columns which we will not need
#fst <- fst[c(1,4:6)]
# Give each site a unique id (for merging with the other datas sets)
#fst$id <- paste(fst$CHROM, fst$pos, sep="_")

# Merge FST data with other results (this will only retain the FST windows in our analysis)
#inconsistent3 <- merge(inconsistent2, fst, by.x=c("chrom", "pos", "id"), by.y=c("CHROM", "pos", "id"), all.x=T)

#summary(inconsistent3$WEIGHTED_FST) # we have 600 windows with no FST values (NA's) - why?
#subset(inconsistent3, is.na(WEIGHTED_FST))$N_VARIANTS # these windows do not have any variants (all NA's) which makes sense why we don't have FST values for them
#sum(subset(inconsistent3, is.na(WEIGHTED_FST))$AIM_count) # 0 - these windows do not have any AIMs (all 0s) which we would expect if there are no variants

## 3. Recombination Rate
#inconsistent_final <- merge(inconsistent3, rcr, by=c("chrom", "pos", "id"), all.x=T)

#inconsistent_unmaskedWallref <- inconsistent_final

# if you did not include AIMs, FST, or rcr for Wall et al. ref
inconsistent_unmaskedWallref <- inconsistent

rm(fst, aims, inconsistent, inconsistent_error, inconsistent2, inconsistent3, calls_record)

# Let's compare the genome-wide distributions of proportion of ancestry state inconsistencies/consistencies for each set of calls
#############################################################################################################################
# Figure S3D
#############################################################################################################################
# inconsistencies
ggplot()  + geom_violin(data=inconsistent_maskedSNPRCref, aes(as.factor(1),prop_inconsistent), fill="mediumpurple3", size=0.75) + geom_violin(data=inconsistent_unmaskedWallref, aes(as.factor(2), prop_inconsistent), fill="indianred3", size=0.75) + geom_boxplot(data=inconsistent_maskedSNPRCref, aes(as.factor(1),prop_inconsistent), fill="white", outlier.color = NA, width=0.5, size=0.75) + geom_boxplot(data=inconsistent_unmaskedWallref, aes(as.factor(2), prop_inconsistent), fill="white", outlier.color = NA, width=0.5, size=0.75) + scale_y_continuous(name="proportion of ancestry state\ninconsistencies across pedigree trios", limits=c(0,1)) + scale_x_discrete(labels=c( "1" = "Masked SNPRC\nreference panel", "2" = "Unmasked Maasai Mara, WNPRC,\nand Mikumi reference panel")) + theme_classic() + theme(text=element_text(size=20), axis.text = element_text(color="black"), axis.title.x = element_blank())

# consistencies (included in the Supplementary Materials)
ggplot()  + geom_violin(data=inconsistent_maskedSNPRCref, aes(as.factor(1),prop_consistent), fill="mediumpurple3", size=0.75) + geom_violin(data=inconsistent_unmaskedWallref, aes(as.factor(2), prop_consistent), fill="indianred3", size=0.75) + geom_boxplot(data=inconsistent_maskedSNPRCref, aes(as.factor(1),prop_consistent), fill="white", outlier.color = NA, width=0.5, size=0.75) + geom_boxplot(data=inconsistent_unmaskedWallref, aes(as.factor(2), prop_consistent), fill="white", outlier.color = NA, width=0.5, size=0.75) + scale_y_continuous(name="proportion of ancestry state\nconsistencies across pedigree trios", limits=c(0,1)) + scale_x_discrete(labels=c( "1" = "Masked SNPRC\nreference panel", "2" = "Unmasked Maasai Mara, WNPRC,\nand Mikumi reference panel")) + theme_classic() + theme(text=element_text(size=20), axis.text = element_text(color="black"), axis.title.x = element_blank()) -> d; d

ggsave("figS3D.png")

#############################################################################################################################
# Model predicting a site's proportion of ancestry state inconsistencies/consistencies across pedigree trios.
#############################################################################################################################

# inconsistencies
summary(lm(prop_inconsistent~AIM_count + WEIGHTED_FST +  log10(n24_anubis), data=inconsistent_maskedSNPRCref))
# consistencies
summary(lm(prop_consistent~AIM_count + WEIGHTED_FST +  log10(n24_anubis), data=inconsistent_maskedSNPRCref))

write.csv(summary(lm(prop_consistent~AIM_count + WEIGHTED_FST + log10(n24_anubis), data=inconsistent_maskedSNPRCref))$coef, "tableS7_results.csv")

#############################################################################################################################
# Proportion of ancestry state inconsistencies/consistencies per pedigree trio.
#############################################################################################################################

# Are some predigree trios more or less inconsistent/consistent than others?
load("local_ancestry_pedigree_trios_maskedSNPRCref.Rd", verbose=T)
load("inconsistent_local_ancestry_calls_pedigree_trios_maskedSNPRCref.Rd", verbose=T) # reload ancestry state inconsistency results generated from the masked, SNPRC reference pabnel

# First, sum the number of inconsistencies (1s) across sites per pedigree trio
tmp <- as.data.frame(colSums(inconsistent[c(3:(nrow(trios)+2))], na.rm=T))
colnames(tmp) <- c("num_inconsistent")
setDT(tmp, keep.rownames = "trio_id")

# Second, sum the number of NAs across sites per pedigree trio
tmp2 <- as.data.frame(colSums(is.na(inconsistent[c(3:(nrow(trios)+2))])))
colnames(tmp2) <- c("num_NA")
setDT(tmp2, keep.rownames = "trio_id")

# Check to make sure everything was calculated correctly
# Check an example trio
tmp3 <- tmp2[10,] # let's pick a random number - 10
tmp3$trio_id -> k
tmp[tmp$trio_id==k,]$num_inconsistent==nrow(subset(inconsistent, AMB_017_AMB_105_AMB_413_1==1)) # TRUE
tmp2[tmp2$trio_id==k,]$num_NA==nrow(subset(inconsistent, is.na(AMB_017_AMB_105_AMB_413_1))) # TRUE

trios_inconsistencies <- merge(tmp, tmp2, by=c("trio_id"))

# Get the number of consistent calls per pedigree trio
trios_inconsistencies$num_consistent <- nrow(inconsistent) - (trios_inconsistencies$num_inconsistent+trios_inconsistencies$num_NA)

# Check to make sure everything was calculated correctly
# Check an example trio
subset(trios_inconsistencies, trio_id=="AMB_315_AMB_111_AMB_135_1")$num_consistent==nrow(subset(inconsistent, AMB_315_AMB_111_AMB_135_1==0)) # TRUE

sum(trios_inconsistencies$num_consistent+trios_inconsistencies$num_inconsistent+trios_inconsistencies$num_NA==nrow(positions))==nrow(trios_inconsistencies) # TRUE

# Get the proportion of inconsistent calls per pedigree trio
trios_inconsistencies$prop_inconsistent <- trios_inconsistencies$num_inconsistent/(trios_inconsistencies$num_inconsistent+trios_inconsistencies$num_consistent)

# Get the proportion of NA calls per pedigree trio
trios_inconsistencies$prop_na <- trios_inconsistencies$num_NA/nrow(inconsistent)

# Get the proportion of consistent calls per pedigree trio
trios_inconsistencies$prop_consistent <- 1-trios_inconsistencies$prop_inconsistent

round(summary(trios_inconsistencies$prop_inconsistent),2)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.06    0.07    0.08    0.09    0.09    0.26 

round(summary(trios_inconsistencies$prop_consistent),2)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.74    0.91    0.92    0.91    0.93    0.94 

# Add coverage and genome-wide anubis ancestry data for each individual in each pedigree trio
# Use the trios data frame which has the data we need
# Make unique trio id to match the trios_inconsistent trio_id in order to merge data frames - should be dad, mom, kid, 1
trios$trio_id <- paste(trios$dad_id, trios$mom_id, trios$kid_id, "1", sep="_")

trios_inconsistencies <- merge(trios_inconsistencies, trios, by=c("trio_id"))
nrow(trios)==nrow(trios_inconsistencies) # TRUE

# Get the minimum coverage across individuals in each pedigree trio
trios_inconsistencies$min_cov <- apply(select(trios_inconsistencies, dad_coverage, mom_coverage, kid_coverage),1,min)

# Check to make sure min_cov was calculated correctly
# Check an example trio
tmp <- trios[29,] # let's pick a random number - 29
tmp$trio_id -> k
min(tmp$kid_coverage, tmp$mom_coverage, tmp$dad_coverage)==trios_inconsistencies[trios_inconsistencies$trio_id==k,]$min_cov # min - TRUE

# inconsistencies
cor.test(trios_inconsistencies$min_cov, trios_inconsistencies$prop_inconsistent)
round(cor.test(trios_inconsistencies$min_cov, trios_inconsistencies$prop_inconsistent)$estimate,3) # -0.553
cor.test(trios_inconsistencies$min_cov, trios_inconsistencies$prop_inconsistent)$p.value #1.067839e-08

# consistencies
cor.test(trios_inconsistencies$min_cov, trios_inconsistencies$prop_consistent)
round(cor.test(trios_inconsistencies$min_cov, trios_inconsistencies$prop_consistent)$estimate,3) # 0.553
cor.test(trios_inconsistencies$min_cov, trios_inconsistencies$prop_consistent)$p.value #1.067839e-08


#############################################################################################################################
# Figure S3E
#############################################################################################################################

# Minimum coverage seems to be the most important factor for predicting ancestry state inconsistencies/consistencies in a pedigree trio
# Plot minimum coverage in a trio (x-axis) vs. proportion of pedigree inconsistencies per trio (y-axis)
ggplot(data=trios_inconsistencies)  + geom_smooth(aes(min_cov, prop_inconsistent), method="lm", color="mediumpurple3") +  scale_x_continuous(name="minimum genome-wide coverage in a pedigree trio", breaks=seq(0,1.4,by=0.2)) + geom_point(aes(min_cov, prop_inconsistent), size=4, alpha=0.7, stroke=1.5, color="grey50") + scale_y_continuous(name="proportion of ancestry state\ninconsistencies in a pedigree trio")  + theme_classic() + theme(text=element_text(size=20), axis.text = element_text(color="black"))

# # Plot minimum coverage in a trio (x-axis) vs. proportion of pedigree consistencies per trio (y-axis) - included in the SI
ggplot(data=trios_inconsistencies)  + geom_smooth(aes(min_cov, prop_consistent), method="lm", color="mediumpurple3") +  scale_x_continuous(name="minimum genome-wide coverage in a pedigree trio", breaks=seq(0,1.4,by=0.2)) + geom_point(aes(min_cov, prop_consistent), size=4, alpha=0.7, stroke=1.5, color="grey50") + scale_y_continuous(name="proportion of ancestry state\nconsistencies in a pedigree trio")  + theme_classic() + theme(text=element_text(size=20), axis.text = element_text(color="black"))

ggsave("figS3E.png")

#save(inconsistent_maskedSNPRCref, inconsistent_unmaskedWallref, trios_inconsistencies, file="inconsistencies_results_latest_20Jul2021.Rd")
#(d | ( e + theme(axis.title.x = element_text(margin = margin(t = -13, unit = "pt")))))
  
