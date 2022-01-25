# Script for recreating figure 1B and S4

# Load R libraries
library(ggplot2)
library(dplyr)
library(data.table)

# Load genotype data needed for PCA 
load("./for_pca.RData", verbose=T) 

#############################################################################################################################
# Figure 1B
#############################################################################################################################

# PCA of high coverage putatively unadmixed individuals and high coverage Amboseli individuals
# Need to remove low coverage data
names$V1[-c(34:46,58:60,65)]
length(names$V1[-c(34:46,58:60,65)])==55 # TRUE - should equal 55

pcrat <- prcomp(covgeno_highcov,scale.=T) # perform PCA on genotype covariance using only high coverage individuals

tmp <- as.data.frame(pcrat$x)

# BGDP anubis individuals are made up of wild caught Aberdares individuals (panu_30877, panu_30977) and SNPRC individuals (panu_L142,  panu_LIV5); BGDP yellow individual is from Mikumi
names_tmp <- names

# find the row numbers of the BGDP individuals you want to relabel
# BGDP individiduals from SNPRC
rownames(names_tmp)[names_tmp$V1 %like% "panu_L"] # 49,50
# BGDP individuals from Aberdares
rownames(names_tmp)[names_tmp$V1 %like% "panu_3"] # 47,48

names_tmp$source <- c(rep("Amboseli",9), rep("SNPRC",24), rep("Mara", 7), rep("WNPRC",6), rep("BGDPanubis",4), rep("SNPRC",7), rep("Mikumi",15))

names_tmp$source[c(49:50)] <- "SNPRC"
names_tmp$source[c(47:48)] <- "Aberdares"

summary(pcrat)

names_tmp$source2 <- factor(names_tmp$source, levels=c("Amboseli", "Mikumi", "Aberdares", "SNPRC")) 

b <- ggplot(data=tmp) + geom_point(aes(PC1, PC2,  col=factor(names_tmp$source2[-c(34:46,58:60,65)]), fill=factor(names_tmp$source2[-c(34:46,58:60,65)])), shape=21, size=4, alpha=0.8, stroke=1.5) + theme_classic() + theme(text=element_text(size=18, family="Helvetica"), legend.position = "none", axis.text = element_text(color="black"))  + scale_fill_manual(values = c("Amboseli"="darkorange2", "Aberdares"="#009E73", "SNPRC"="grey70", "Mikumi"="gold1")) + scale_color_manual(values = c("Amboseli"="darkorange2", "Aberdares"="springgreen4", "SNPRC"="grey55", "Mikumi"="gold2")) + scale_x_continuous(name="PC1 (84% variance explained)") + scale_y_continuous(name="PC2 (2% variance explained)"); b
# get legend
ggplot(data=tmp) + geom_point(aes(PC1, PC2,  col=factor(names_tmp$source2[-c(34:46,58:60,65)]), fill=factor(names_tmp$source2[-c(34:46,58:60,65)])), shape=21, size=4, alpha=0.8, stroke=1.5) + theme_classic() + theme(text=element_text(size=25), legend.position = "right", legend.title = element_blank(), axis.text = element_text(color="black"))  + scale_fill_manual(values = c("Amboseli"="darkorange2", "Aberdares"="#009E73", "SNPRC"="grey70", "Mikumi"="gold1")) + scale_color_manual(values = c("Amboseli"="darkorange2", "Aberdares"="springgreen4", "SNPRC"="grey55", "Mikumi"="gold2")) + scale_x_continuous(name="PC1 (0.84)") + scale_y_continuous(name="PC2 (0.02)")

ggsave("fig1B.png")

# zoom in on "yellow-like" individuals
ggplot(data=tmp) + geom_point(aes(factor(as.factor(names_tmp$source2[-c(34:46,58:60,65)]), levels=c("Aberdares", "Amboseli","SNPRC","Mikumi")), PC1, col=factor(names_tmp$source2[-c(34:46,58:60,65)]), fill=factor(names_tmp$source2[-c(34:46,58:60,65)])), size=4, alpha=0.8, shape=21, stroke=1.5) + coord_flip() + theme_classic() + theme(text=element_text(size=28, family="Helvetica"), legend.position = "none", axis.text = element_text(color="black"), axis.title.y = element_blank(), axis.text.y = element_text(angle=30, vjust=-1), axis.ticks.y = element_blank()) + scale_fill_manual(values = c("Amboseli"="darkorange2", "Aberdares"="#009E73", "SNPRC"="grey70", "Mikumi"="#FFED4F")) + scale_color_manual(values = c("Amboseli"="darkorange2", "Aberdares"="springgreen4", "SNPRC"="grey55", "Mikumi"="gold1")) + scale_y_continuous(limits=c(-7.8,-6)) + scale_x_discrete(limits=c("Amboseli", "SNPRC", "Mikumi"))

ggsave("fig1B_inset.png")

# Combine 1B and 1C for Figure 1
#library(patchwork)
#(b + theme(plot.margin = unit(c(0,20,0,0), "pt"))) + (c + theme(plot.margin = unit(c(0,0,0,20), "pt")))

#############################################################################################################################
# Figure S4
#############################################################################################################################

# PCA of all putatively unadmixed individuals (both low and high coverage) and high coverage Amboseli individuals
pcrat <- prcomp(covgeno,scale.=T) # perform PCA on genotype covariance using all individuals regardless of their coverage

tmp <- as.data.frame(pcrat$x)

# BGDP anubis individuals are made up of wild caught Aberdares individuals (panu_30877, panu_30977) and SNPRC individuals (panu_L142,  panu_LIV5); BGDP yellow individual is from Mikumi
names_tmp <- names
names_tmp$source <- c(rep("Amboseli",9), rep("SNPRC",24), rep("Mara", 7), rep("WNPRC",6), rep("BGDPanubis",4), rep("SNPRC",7), rep("Mikumi",15))

# find the row numbers of the BGDP individuals you want to relabel
# BGDP individiduals from SNPRC
rownames(names_tmp)[names_tmp$V1 %like% "panu_L"] # 49,50
# BGDP individuals from Mikumi
rownames(names_tmp)[names_tmp$V1 %like% "pcyn_16066"] # 72
# BGDP individuals from Aberdares
rownames(names_tmp)[names_tmp$V1 %like% "panu_3"] # 47,48

names_tmp$source[c(49:50)] <- "SNPRC"
names_tmp$source[c(47:48)] <- "Aberdares"

summary(pcrat)

names_tmp$source2 <- factor(names_tmp$source, levels=c("Amboseli", "Mikumi", "Aberdares", "Mara", "WNPRC", "SNPRC")) 
# assign low and high coverage
names_tmp$coverage <- NA
names_tmp$coverage[-c(34:46,58:60,65)] <- "high"
names_tmp[,c("coverage")][is.na(names_tmp[,c("coverage")])] <- "low"

# without legend
ggplot(data=tmp) + geom_point(aes(PC1, PC2,  col=factor(names_tmp$source2), fill=factor(names_tmp$source2), shape=factor(names_tmp$coverage)), size=4, alpha=0.8, stroke=1.5) + theme_classic() + theme(text=element_text(size=18), axis.text = element_text(color="black"), legend.position = "none")  + scale_fill_manual(values = c("Amboseli"="darkorange2", "Aberdares"="#009E73", "SNPRC"="grey70", "Mikumi"="gold1", "WNPRC"="forestgreen", "Mara"="#00B81F")) + scale_color_manual(values = c("Amboseli"="darkorange2", "Aberdares"="#009E73", "SNPRC"="grey70", "Mikumi"="gold1", "WNPRC"="forestgreen", "Mara"="#00B81F"))  + scale_shape_manual(values=c("high"=21, "low"=24)) + scale_x_continuous(name="PC1 (0.53 PVE)") + scale_y_continuous(name="PC2 (0.16 PVE)")
# with legend
ggplot(data=tmp) + geom_point(aes(PC1, PC2,  col=factor(names_tmp$source2), fill=factor(names_tmp$source2), shape=factor(names_tmp$coverage)), size=4, alpha=0.8, stroke=1.5) + theme_classic() + theme(text=element_text(size=18), axis.text = element_text(color="black"))  + scale_fill_manual(values = c("Amboseli"="darkorange2", "Aberdares"="#009E73", "SNPRC"="grey70", "Mikumi"="gold1", "WNPRC"="forestgreen", "Mara"="#00B81F")) + scale_color_manual(values = c("Amboseli"="darkorange2", "Aberdares"="#009E73", "SNPRC"="grey70", "Mikumi"="gold1", "WNPRC"="forestgreen", "Mara"="#00B81F"))  + scale_shape_manual(values=c("high"=21, "low"=24)) + scale_x_continuous(name="PC1 (0.53 PVE)") + scale_y_continuous(name="PC2 (0.16 PVE)")

ggsave("figS4.png")
