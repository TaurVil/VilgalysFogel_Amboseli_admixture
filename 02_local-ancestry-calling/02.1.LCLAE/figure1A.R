# Script for recreating figure 1A local ancestry tracts

# Load R libraries
library(ggplot2)
library(dplyr)

# Plot local ancestry estimates for multiple individuals where a horizontal row corresponds to the first 20 Mb of chromosome 1 for each individual, organized from top to bottom from the most anubis to most yellow individuals

# Below is code for plotting local ancestry tracts for Amboseli individuals but can be expanded to the other populations included in figure 1A by replacing the tracts file below with the tracts files from those populations (figure 1A plots ancestry estimates using unmasked data for non-Amboseli individuals)
# Load ancestry tracts for Amboseli individuals
tracts <- read.table("amboseli_LCLAE_tracts.txt", header=T)

# Remove tracts <1 kb in length
tracts <- tracts[tracts$length2>1000,]

# Change ancestry state into a factor variable
tracts$state2 <- as.factor(tracts$state)

# Calculate genome-wide anubis ancestry so we can organize individuals from the most anubis to the most yellow (this produces the same genome-wide anubis ancestry estimates as in Table S1)
overall_admix <- tracts %>% group_by(table_s1_id, state) %>% summarize(total_state_length=sum(length2), .groups="keep") %>% as.data.frame()
overall_admix2 <- overall_admix %>% group_by(table_s1_id) %>% mutate(total_length=sum(total_state_length)) %>% as.data.frame()
overall_admix2$state_length_multiple <- overall_admix2$state*overall_admix2$total_state_length
overall_admix3 <- overall_admix2 %>% group_by(table_s1_id) %>% mutate(anubis_admix=sum(state_length_multiple)/(2*total_length)) %>% as.data.frame()

overall_admix4 <- overall_admix3[c(1,4,6)]
overall_admix5 <- unique(overall_admix4[c("table_s1_id", "total_length", "anubis_admix")])

tracts2 <- merge(tracts, overall_admix5, by.x=c("table_s1_id"), by.y=c("table_s1_id"))
tracts2$table_s1_id_reordered <- reorder(tracts2$table_s1_id, tracts2$anubis_admix)

overall_admix_amb <- overall_admix5
tracts_amb <- tracts2

# For Amboseli, plot a random subsample of 100 individuals
set.seed(1234) # so the random subsampling is repeatale
random_indiv <- sample(unique(tracts_amb$table_s1_id_reordered), 100, replace=F)
tracts_amb_random <- tracts_amb[tracts_amb$table_s1_id_reordered %in% random_indiv,]
tracts_amb_random$table_s1_id_reordered <- droplevels(tracts_amb_random$table_s1_id_reordered)

# Visualize a subset of the genome
# Let's look at the first 20 Mb of chromosome 1. Note that most tracts do not switch ancestry states exactly at the 20 Mb position of chromosome 1 so subset the data and visualize tracts that extend beyond 20 Mb (let's go to 30 Mb) but draw a 20 Mb line so we can see the end of the genomic region we're interested in
ggplot(subset(tracts_amb_random, chrom=="chr1" & end<=30000000)) + geom_rect(aes(xmin = as.numeric(table_s1_id_reordered) - 0.2, xmax = as.numeric(table_s1_id_reordered) + 0.2, ymin = start, ymax = end, color=state2, fill=state2)) + coord_flip() + scale_x_continuous(name="individual",limits=c(0,101), breaks=seq(0,101,by=50)) + scale_y_continuous(labels = scales::comma) + scale_color_manual(values=c("#FFED4F", "orange2", "#009E73"), name="ancestry assignment",labels=c("homozygous yellow", "heterozygous yellow-anubis", "homozygous anubis")) + scale_fill_manual(values=c("#FFED4F", "orange2", "#009E73"), name="ancestry assignment",labels=c("homozygous yellow", "heterozygous yellow-anubis", "homozygous anubis")) + theme_void() + theme(text=element_text(size=25), legend.title = element_blank(), legend.position = "none") + geom_hline(aes(yintercept=20000000)) 

ggsave("fig1A_amboseli.png")

rm(overall_admix, overall_admix2, overall_admix3, overall_admix4, overall_admix5, overall_admix, tracts, tracts2, random_indiv)
