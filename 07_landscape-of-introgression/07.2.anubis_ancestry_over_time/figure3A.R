# Script for recreating figure 3A

# Load R libraries
library(ggplot2)
library(dplyr)

# Load demographic data for Amboseli individuals with confirmed ids, including the starting and endings years they were present in the population
demog <- read.table("amboseli.demographic.info.txt", header=T)

# Get a list of individuals in the population that were present in each year
years_list <- seq(min(demog$start_year), max(demog$end_year)) # get list of all unique years in the demographic data
max(years_list) # 2021
# Remove 2021 because we do not have data for the entire year
years_list <- years_list[1:length(years_list)-1]

# Store results of mean genome-wide anubis ancestry in the Amboseli baboon population each year
years <- data.frame()

for (i in 1:length(years_list)) { # do for each year
  years[i,1] <- years_list[i] 
  
  # Get all individuals to be included in the calculation for the focal year
  tmp <- demog[demog$start_year<=years_list[i] & demog$end_year>=years_list[i],]
  years[i,2] <- mean(tmp$genome_wide_anubis_ancestry) # get the mean anubis ancestry in the focal year for all individuals in the population
  years[i,3] <- nrow(tmp) # get the sample size used in the calculation (should also be the same number of rows as in tmp)
  
  # Get all immigrants (a subset of all individuals) to be included in the calculation for the focal year
  tmp2 <- tmp[tmp$entry_type=="I",]
  years[i,4] <- mean(tmp2$genome_wide_anubis_ancestry) # get the mean anubis ancestry in the focal year for all immigrants in the population
  years[i,5] <- nrow(tmp[tmp$entry_type=="I",]) # get the sample size used in the calculation (should also be the same number of rows as in tmp2)
  
  # Get all individuals who are not immigrants (i.e., they were born into study groups or present at the beginning of observations; again, a subset of all individuals) to be included in the calculation for the focal year
  tmp3 <- tmp[tmp$entry_type!="I",]
  years[i,6] <- mean(tmp3$genome_wide_anubis_ancestry) # get the mean anubis ancestry in the focal year for all non-immigrants in the population
  years[i,7] <- nrow(tmp[tmp$entry_type!="I",]) # get the sample size used in the calculation (should also be the same number of rows as in tmp3)
  
  # Get all males (a subset of all individuals) to be included in the calculation for the focal year
  tmp4 <- tmp[tmp$sex=="M",] 
  years[i,8] <- mean(tmp4$genome_wide_anubis_ancestry) # get the mean anubis ancestry in the focal year for all males in the population
  years[i,9] <- nrow(tmp[tmp$sex=="M",]) # get the sample size used in the calculation (should also be the same number of rows as in tmp4)
  
  # Get all females (a subset of all individuals) to be included in the calculation for the focal year
  tmp5 <- tmp[tmp$sex=="F",] 
  years[i,10] <- mean(tmp5$genome_wide_anubis_ancestry) # get the mean anubis ancestry in the focal year for all females in the population
  years[i,11] <- nrow(tmp[tmp$sex=="F",]) # get the sample size used in the calculation (should also be the same number of rows as in tmp5)
  
  print(i)
}

colnames(years) <- c("year", "mean_anubis_ancestry", "n_total", "mean_anubis_ancestry_immigrant", "n_immigrant", "mean_anubis_ancestry_notimmigrant", "n_not_immigrant", "mean_anubis_ancestry_male", "n_male", "mean_anubis_ancestry_female", "n_female")

# Plot the mean anubis ancestry per year for all individuals, immigrants, and non-immigrants
ggplot(data=years) + geom_line(aes(year, mean_anubis_ancestry), size=3, alpha=0.8) + geom_line(aes(year, mean_anubis_ancestry_immigrant), size=3, color="#009E73", alpha=0.8) + geom_line(aes(year, mean_anubis_ancestry_notimmigrant), size=3, color="grey50", alpha=0.8) + geom_text(aes(year, y=0.283, label=n_total),size=4, angle=45, fontface=c("bold")) + geom_text(aes(year, y=0.27, label=n_immigrant),size=4,angle=45, color="#009E73", fontface=c("bold")) + scale_y_continuous(name="mean genome-wide anubis ancestry")  + theme_classic() + theme(text=element_text(size=20, family="Helvetica"), axis.text = element_text(color="black")) + geom_point(aes(2004, y=0.259), color="#88A055", pch=8, size=3, stroke=1.5) +  geom_point(aes(2001, y=0.259), color="#FACF38", pch=8, size=3, stroke=1.5)  + geom_point(aes(2006, y=0.259), color="#019E73", pch=8, size=3, stroke=1.5) 

ggsave("fig3A_main.pdf")

# Plot the mean anubis ancestry per year for females and males
ggplot(data=years) + geom_line(aes(year, mean_anubis_ancestry_male), size=3, alpha=0.8, color="#0A81FF") + geom_line(aes(year, mean_anubis_ancestry_female), size=3, alpha=0.8, color="#B82862") +  scale_y_continuous(name="mean genome-wide\nanubis ancestry", limits=c(0.27, 0.5), breaks=seq(0.3,0.5, by=0.05)) + theme_classic() + theme(text=element_text(size=38, family="Helvetica"), axis.text = element_text(color="black"))

ggsave("fig3A_inset.pdf")
