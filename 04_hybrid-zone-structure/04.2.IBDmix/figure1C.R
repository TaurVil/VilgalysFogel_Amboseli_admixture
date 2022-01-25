# Script for recreating figure 1C

# Load R libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Load IBDmix results for putative unadmixed yellow individuals (Mikumi, SNPRC yellow). First column includes the individual ID and second column includes the individual's population. The remaining columns represents the mean estimated IBD between the individual (row) and the population (column name). 
mikumi_SNPRCyellow_ibd <- read.table("ibdmix_yellow_estimates.txt", header=T) 
# Load IBDmix results for high coverage Amboseli individuals
amboseli_ibd <- read.table("ibdmix_amboseli_estimates.txt", header=T) 

#############################################################################################################################
# Figure 1C
#############################################################################################################################

# Plot the IBDmix results using boxplots and jittered points
set.seed(1234) # for reproducibility of jittered points

# Color anubis, hamadryas, guinea following the color scheme in Fig. 1A map (green = anubis, blue = hamadryas, purple = guinea)
c <- ggplot() + geom_jitter(data=mikumi_SNPRCyellow_ibd[mikumi_SNPRCyellow_ibd$population=="Mikumi",], aes(x=1, y=SNPRCanubis), alpha=0.8, width=0.5, color="#009E73", fill="#009E73", size=4, shape=21, stroke=1.5) + geom_boxplot(data=mikumi_SNPRCyellow_ibd[mikumi_SNPRCyellow_ibd$population=="Mikumi",], aes(x=1, y=SNPRCanubis), alpha=0.7, outlier.color=NA, fill="#009E73", size=0.75) + geom_jitter(data=mikumi_SNPRCyellow_ibd[mikumi_SNPRCyellow_ibd$population=="Mikumi",], aes(x=2, y=hamadryas), alpha=0.8, width=0.5, color="dodgerblue3", fill="dodgerblue3",  size=4, shape=21, stroke=1.5) + geom_boxplot(data=mikumi_SNPRCyellow_ibd[mikumi_SNPRCyellow_ibd$population=="Mikumi",], aes(x=2, y=hamadryas), alpha=0.7, outlier.color=NA, fill="dodgerblue3", size=0.75)  +  geom_jitter(data=mikumi_SNPRCyellow_ibd[mikumi_SNPRCyellow_ibd$population=="Mikumi",], aes(x=3, y=guinea), alpha=0.8, width=0.5, color="mediumorchid3",fill="mediumorchid3",  size=4, shape=21, stroke=1.5) + geom_boxplot(data=mikumi_SNPRCyellow_ibd[mikumi_SNPRCyellow_ibd$population=="Mikumi",], aes(x=3, y=guinea), alpha=0.7, outlier.color=NA, fill="mediumorchid3", size=0.75) + geom_jitter(data=mikumi_SNPRCyellow_ibd[mikumi_SNPRCyellow_ibd$population=="SNPRCyellow",], aes(x=7, y=SNPRCanubis), alpha=0.8, width=0.5, color="#009E73",fill="#009E73",  size=4, shape=21, stroke=1.5) + geom_boxplot(data=mikumi_SNPRCyellow_ibd[mikumi_SNPRCyellow_ibd$population=="SNPRCyellow",], aes(x=7, y=SNPRCanubis), alpha=0.7, outlier.color=NA, fill="#009E73", size=0.75)  +  geom_jitter(data=mikumi_SNPRCyellow_ibd[mikumi_SNPRCyellow_ibd$population=="SNPRCyellow",], aes(x=8, y=hamadryas), alpha=0.8, width=0.5, color="dodgerblue3", fill="dodgerblue3", size=4, shape=21, stroke=1.5) + geom_boxplot(data=mikumi_SNPRCyellow_ibd[mikumi_SNPRCyellow_ibd$population=="SNPRCyellow",], aes(x=8, y=hamadryas), alpha=0.7, outlier.color=NA, fill="dodgerblue3", size=0.75) + geom_jitter(data=mikumi_SNPRCyellow_ibd[mikumi_SNPRCyellow_ibd$population=="SNPRCyellow",], aes(x=9, y=guinea), alpha=0.8, width=0.5, color="mediumorchid3", fill="mediumorchid3",  size=4, shape=21, stroke=1.5)  + geom_boxplot(data=mikumi_SNPRCyellow_ibd[mikumi_SNPRCyellow_ibd$population=="SNPRCyellow",], aes(x=9, y=guinea), alpha=0.7, outlier.color=NA, fill="mediumorchid3", size=0.75)  +
  geom_jitter(data=amboseli_ibd, aes(x=13, y=SNPRCanubis), alpha=0.8, width=0.5, color="#009E73",  fill="#009E73", size=4, shape=21, stroke=1.5) + geom_boxplot(data=amboseli_ibd, aes(x=13, y=SNPRCanubis), alpha=0.7, outlier.color=NA, fill="#009E73", size=0.75) +  geom_jitter(data=amboseli_ibd, aes(x=14, y=hamadryas), alpha=0.8, width=0.5, color="dodgerblue3", fill="dodgerblue3", size=4, shape=21, stroke=1.5) + geom_boxplot(data=amboseli_ibd, aes(x=14, y=hamadryas), alpha=0.7, outlier.color=NA, fill="dodgerblue3", size=0.75) + geom_jitter(data=amboseli_ibd, aes(x=15, y=guinea), alpha=0.8, width=0.5, color="mediumorchid3", fill="mediumorchid3", size=4, shape=21, stroke=1.5)  + geom_boxplot(data=amboseli_ibd, aes(x=15, y=guinea), alpha=0.7, outlier.color=NA, fill="mediumorchid3", size=0.75)  +  theme_classic() + theme(text=element_text(size=18, family="Helvetica"), axis.text = element_text(color="black")) + scale_x_continuous("population", breaks=c(3,8,14), labels=c("Mikumi", "SNPRC yellow founders", "Amboseli")) + scale_y_continuous("proportion identical-by-descent"); c

# For the legend
# need to reshape the data for easier plotting
ibd2 <- mikumi_SNPRCyellow_ibd %>% pivot_longer(c("SNPRCanubis", "guinea", "hamadryas", "aberdares"), names_to = "species", values_to = "ancestry") 
set.seed(1234)
ggplot(data=ibd2[ibd2$species!="aberdares",])  + geom_jitter(aes(x=population, y=ancestry, color=species), alpha=0.8, width=0.05, size=2, stroke=1) + scale_fill_manual(values=c("#009E73","dodgerblue3",  "mediumorchid3"), labels=c("anubis", "hamadryas", "guinea")) + scale_color_manual(values=c("#009E73","dodgerblue3", "mediumorchid3"), labels=c("anubis",  "hamadryas", "guinea")) + theme_classic() + theme(text=element_text(size=25), axis.text = element_text(color="black")) + scale_y_continuous("proportion identical-by-descent")

ggsave("fig1C.png")
