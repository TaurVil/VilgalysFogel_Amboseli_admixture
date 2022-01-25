# Script for recreating figure 3C

# Load R libraries
library(ggplot2)
library(dplyr)
library(scales)

# Plot genome-wide anubis ancestry for all Amboseli individuals (n=442)
# Load Table S1 data
load("VilgalysFogel_main_data_file.250kb_windows.RData")

# Plot histogram with density lines for historic vs. recent hybrids
ggplot() + geom_histogram(data=ids[ids$population_or_species=="Amboseli",], aes(x=genome_wide_anubis_ancestryd, fill=..x..), color="white", binwidth=0.005) + geom_density(data=ids[ids$population_or_species=="Amboseli" & ids$recent_hybridsc=="historical",], aes(x=genome_wide_anubis_ancestryd), color="grey60", size=1.5) +  geom_density(data=ids[ids$population_or_species=="Amboseli" & ids$recent_hybridsc=="recent",], aes(x=genome_wide_anubis_ancestryd), color="grey60", size=1.5) + scale_x_continuous(name="genome-wide anubis ancestry") + scale_y_continuous(name="count") + theme_classic() + theme(text=element_text(size=30, family="Helvetica"), legend.position = "none", axis.title.x = element_text(vjust=-0.02), plot.margin = unit(c(0,0,10,0), "pt")) + scale_fill_gradientn(colors=c("#FFED4F", "orange","#009E73"), values=scales::rescale(c(0,0.25,0.4,0.5,1)))+ coord_cartesian(xlim=c(0,1))
ggsave("fig3C.png")
