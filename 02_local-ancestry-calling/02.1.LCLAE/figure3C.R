# Script for recreating figure 3C

# Load R libraries
library(ggplot2)
library(dplyr)
library(scales)

# Plot genome-wide anubis ancestry for all Amboseli individuals (n=442)
# Load Table S1 data
sample_info <- read.table("XXXXXX") # Tauras, I'm Assuming Table S1 will be in the final Rdata

# Plot histogram
ggplot(sample_info, aes(x=genome_wide_anubis_ancestryd, fill=..x..)) + geom_histogram(color="white", binwidth=0.005) + scale_x_continuous(name="genome-wide anubis ancestry") + theme_classic() + theme(text=element_text(size=25), legend.position = "none", axis.title.x = element_text(vjust=-0.02), plot.margin = unit(c(0,0,10,0), "pt")) + scale_fill_gradientn(colors=c("#FFED4F", "orange","#009E73"), values=scales::rescale(c(0,0.25,0.4,0.5,1)))+ coord_cartesian(xlim=c(0,1))

ggsave("fig3C.png")
