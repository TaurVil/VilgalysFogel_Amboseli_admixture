# Script for recreating figure S6 "Estimated timing of admixture in Amboseli baboons"

# Load R libraries
library(ggplot2)

# Load DATES results
dates <- read.table("DATES_results.txt", header=T)

# Load ancestry tracts generated using masked SNPRC individuals for the reference panel (used in main text analyses) available in the Duke Data Repository at XXX
tracts <- read.table("amboseli.tracts.maskedSNPRCref.txt", header=T)

# Remove tracts <1 kb in length
tracts <- tracts[tracts$length2>1000,]

#############################################################################################################################
# Figure 6A
#############################################################################################################################

# Plot DATES estimates for recent (blue/purple) vs. historical (brown) hybrids with high coverage data
ggplot(data=dates, aes(x=hybrid_status, y=DATES_estimate)) +
    geom_boxplot(alpha=0.5, outlier.color=NA,  size=0.75,  aes(fill=hybrid_status)) + 
    geom_point(alpha=0.95, size=4, shape=21, stroke=1.5, aes(fill=hybrid_status)) +
    theme_classic() +
    theme(text=element_text(size=18), axis.text = element_text(color="black"), legend.position = "none") +
    scale_x_discrete(name="hybrid status") +
    scale_y_continuous("generations since admixture") +
    scale_fill_manual(values=c("#AD7748", "#5B507A"))

ggsave("figS6A.png")

#############################################################################################################################
# Figure 6B
#############################################################################################################################

# Plot the tract lengths for two representative individuals in panel A
mean(dates[dates$hybrid_status=="historical",]$DATES_estimate) # for historically admixed individuals, let's use the mean = 283.1643 (which is very close to AMB_310's estimate of 316.12
# for recently admixed individuals, let's use AMB_301	

# define our tract types for plotting
tract_types <- c(`2` = "homozygous anubis\n ancestry tracts", `1`="heterozygous yellow-anubis\nancestry tracts", `0`="homozygous yellow\n ancestry tracts")

# Get tract lengths for the the historically and recently admixed examples individuals
tmp <- tracts[tracts$table_s1_id=="AMB_310" | tracts$table_s1_id=="AMB_301",]

tmp$table_s1_id <-droplevels(tmp$table_s1_id)

# Add a factor variable version of ancestry state for plotting
tmp$state2 <- as.factor(tmp$state)

# Plot log 10 of tract lengths colored by the hybrid status (recent (blue/purple) vs. historical (brown)) of the individual
ggplot(tmp, aes(x=log10(length2), fill=table_s1_id, alpha=0.90)) + geom_density() + facet_grid(as.factor(state2)~., labeller = as_labeller(tract_types)) + scale_x_continuous(name="log10(ancestry tract length)") +  theme_classic() + theme(axis.text=element_text(size=25, colour = "black"),axis.title = element_text(size=25), strip.text=element_text(size=25), legend.position = "none") + scale_fill_manual(values=c("#5B507A", "#AD7748")) 

ggsave("figS6B.png")
