# Script for recreating figure 3B ancestry tracts
# Pedigree schematic and additional changes to figure (e.g., formatting of baboon female photos, asterisks, etc.) were done in keynote 8.1

# Load R libraries
library(ggplot2)

# Load ancestry tracts for Amboseli individuals
tracts <- read.table("amboseli_LCLAE_tracts.txt", header=T)

# Remove tracts <1 kb in length
tracts <- tracts[tracts$length2>1000,]

# Change ancestry state into a factor variable
tracts$state2 <- as.factor(tracts$state)

# Order the chromosome variable
tracts$chrom_reordered <- factor(tracts$chrom, levels = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20"))

# We'll plot ancestry tracts genome-wide for example Amboseli females of historic (AMB_202) and recent (AMB_044) hybrid ancestry. 
# Individual with historic ancestry (AMB_202)
left <- ggplot(subset(tracts, table_s1_id=="AMB_202")) + geom_rect(aes(xmin = as.numeric(table_s1_id) - 1, xmax = as.numeric(table_s1_id) + 1, ymin = start, ymax = end, color=state2, fill=state2)) + coord_flip() + facet_grid(chrom_reordered~.) +  scale_y_continuous(labels = scales::comma) + scale_color_manual(values=c("#FFED4F", "orange2", "#009E73"), name="ancestry assignment",labels=c("homozygous yellow", "heterozygous yellow-anubis", "homozygous anubis")) + scale_fill_manual(values=c("#FFED4F", "orange2", "#009E73"), name="ancestry assignment",labels=c("homozygous yellow", "heterozygous yellow-anubis", "homozygous anubis")) + theme_classic() + theme(text=element_text(size=25), axis.ticks = element_blank(), legend.position = "none", axis.text = element_text(color="black"), axis.text.y=element_blank()); left

# Individual with recent ancestry (AMB_044)
right <- ggplot(subset(tracts, table_s1_id=="AMB_044")) + geom_rect(aes(xmin = as.numeric(table_s1_id) - 1, xmax = as.numeric(table_s1_id) + 1, ymin = start, ymax = end, color=state2, fill=state2)) + coord_flip() + facet_grid(chrom_reordered~.) +  scale_y_continuous(labels = scales::comma) + scale_color_manual(values=c("#FFED4F", "orange2", "#009E73"), name="ancestry assignment",labels=c("homozygous yellow", "heterozygous yellow-anubis", "homozygous anubis")) + scale_fill_manual(values=c("#FFED4F", "orange2", "#009E73"), name="ancestry assignment",labels=c("homozygous yellow", "heterozygous yellow-anubis", "homozygous anubis")) + theme_classic() + theme(text=element_text(size=25), axis.ticks = element_blank(), legend.position = "none", axis.text = element_text(color="black"), axis.text.y=element_blank()); right

# Plot both individuals
left | right

ggsave("fig3B.png") # cropped out x and y axis labels later in keynote
