# Script for recreating figures 2 and 3D-F

# Load R libraries
library(ggplot2)
library(gtools)
library(ggbeeswarm)
library(patchwork)

#############################################################################################################################
# Figure 2
#############################################################################################################################

# Load baboon results
load("../VilgalysFogel_main_data_file.250kb_windows.RData", verbose=T)

# Rename "to_analyze" data frame to "baboon"
baboon <- to_analyze

# For plotting purposes only, assign fixed differences, B values, and recombination rate to quantiles
baboon$fixed_quant <- quantcut(baboon$fixed, q=5) # assign fixed differences to quantiles
baboon$B_quant <- quantcut(baboon$B, q=5) # assign B values to quantiles
baboon$rcr_quant <- quantcut(baboon$rcr, q=5) # assign recombination rate to quantiles

#############################################################################################################################
# Figure 2A, 2C, 2E
#############################################################################################################################

# Get the median anubis ancestry across all 250 kb windows so we can plot a dashed line in each panel
median_anubis_ancestry <- median(baboon$mean_ancestry)
round(median_anubis_ancestry,3) # 0.367

# Use swarm plots and boxplots for visualizing patterns in the baboon data
# Plot fixed differences vs. introgressed anubis ancestry
plotA <- ggplot(data = baboon, aes(x=fixed_quant, y=mean_ancestry)) + 
  geom_beeswarm(aes(color=as.numeric(fixed_quant)), cex=0.4, alpha=0.5) + 
  geom_boxplot(width=0.15,outlier.shape = NA, size=1) +
  geom_hline(aes(yintercept=median_anubis_ancestry), linetype="dashed", size=0.75, color="grey50") +
  theme_classic() + 
  theme(text=element_text(size=14), axis.text = element_text(color="black"), legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_text(vjust=-2), plot.margin = margin(b=10))  + 
  labs(x="fixed differences quantile", y = "") + scale_color_gradient(low="#B09FD1", high="#5C438A"); plotA

# Plot B values vs. introgressed anubis ancestry
plotC <- ggplot(data = baboon, aes(x=B_quant, y=mean_ancestry)) + 
  geom_beeswarm(aes(color=as.numeric(B_quant)), cex=0.4, alpha=0.5) + 
  geom_boxplot(width=0.15,outlier.shape = NA, size=1) + 
  geom_hline(aes(yintercept=median_anubis_ancestry), linetype="dashed", size=0.75, color="grey50") +
  theme_classic() + 
  theme(text=element_text(size=14), axis.text = element_text(color="black"), legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_text(vjust=-2), plot.margin = margin(b=10))  +
  labs(x="B value quantile", y = "proportion of introgressed (anubis) ancestry") + scale_color_gradient(low="#9CAED3", high="#506FB3"); plotC

# Plot recombination rate vs. introgressed anubis ancestry
plotE <- ggplot(data = baboon, aes(x=rcr_quant, y=mean_ancestry)) + 
  geom_beeswarm(aes(color=as.numeric(rcr_quant)), cex=0.4, alpha=0.5)  + 
  geom_boxplot(width=0.15,outlier.shape = NA, size=1) +
  geom_hline(aes(yintercept=median_anubis_ancestry), linetype="dashed", size=0.75, color="grey50") +
  theme_classic() + 
  theme(text=element_text(size=14), axis.text = element_text(color="black"), legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_text(vjust=-2), plot.margin = margin(b=10)) + 
  labs(x="recombination rate quantile", y = "") + scale_color_gradient(low="#9BC6D4", high="#366B7D"); plotE


#############################################################################################################################
# Figure 2B, 2D, 2F
#############################################################################################################################

# Load hominin data
load("./windows.human.250kb.RData", verbose=T)
rm(distance, Neand_ancestry)

# Rename "all" data frame to "hominin"
hominin <- all

# Only retain windows where more than 60% of the 250 kb window includes "testable" parts of the human genome (i.e., they are not regions labeled as difficult to detect Neanderthal ancestry and in regions where Neanderthal genotypes have been called)
hominin <- hominin[hominin$testable > 150000,] 
# Only retain windows where we have non-NA values for all data
hominin <- hominin[rowSums(is.na(hominin)) == 0,] # we lose one line

# As with the baboon results above, for plotting purposes only, assign fixed differences, B values, and recombination rate to quantiles
hominin$fixed_N3_quant <- quantcut(hominin$fixed_N3, q=5) # assign fixed differences to quantiles
hominin$B_quant <- quantcut(hominin$B, q=5) # assign b-values to quantiles
hominin$rcr_quant <- quantcut(hominin$Icelandic, q=5) # assign recombination rate to quantiles

# To compare slopes between species (baboons and hominins) with different levels of introgressed ancestry, center mean introgressed ancestry on 0 within species and divide by the sd to facilitate comparison of slopes
# Mean center the baboon data
baboon$mean_ancestry_centered <- baboon$mean_ancestry - mean(baboon$mean_ancestry)
# Mean center the hominin data
hominin$mean_ancestry_centered <- hominin$mean_ancestry - mean(hominin$mean_ancestry)

# Standardize the mean centered data by the sd in ancestry
# Divide the mean centered baboon data by the sd
baboon$mean_ancestry_centered_sd <- baboon$mean_ancestry_centered/ sd(baboon$mean_ancestry)
# Divde the mean centered hominin data by the sd
hominin$mean_ancestry_centered_sd <- hominin$mean_ancestry_centered/ sd(hominin$mean_ancestry)

# Plot fixed differences vs. introgressed ancestry for hominins and baboons
plotB <- ggplot() + 
  geom_smooth(data = hominin, aes(x=rank(fixed_N3)/length(fixed_N3), y=mean_ancestry_centered_sd), color="#5C438A", size=1.5, linetype = "longdash", method="lm", se=FALSE) + geom_smooth(data = baboon, aes(x=rank(fixed)/length(fixed), y=mean_ancestry_centered_sd), color="#5C438A", size=1.5, method="lm", se=FALSE) + theme_classic() +
  theme(text=element_text(size=14), axis.text = element_text(color="black"), axis.text.x=element_blank(),axis.ticks.x=element_blank(), legend.position = "none", axis.title.x = element_text(vjust=-2), plot.margin = margin(l=20, b=10))  + 
  labs(x="rank-ordered number of fixed differences", y=""); plotB

# Plot B values vs. introgressed ancestry for hominins and baboons
plotD <- ggplot() + 
  geom_smooth(data = hominin, aes(x=rank(B)/length(B), y=mean_ancestry_centered_sd), color="#506FB3", size=1.5, linetype = "longdash", method="lm", se=FALSE) + geom_smooth(data = baboon, aes(x=rank(B)/length(B), y=mean_ancestry_centered_sd), color="#506FB3", size=1.5, method="lm", se=FALSE) + theme_classic() +
  theme(text=element_text(size=14), axis.text = element_text(color="black"), axis.text.x=element_blank(),axis.ticks.x=element_blank(), legend.position = "none", axis.title.x = element_text(vjust=-2), plot.margin = margin(l=20, b=10))  + 
  labs(x="rank-ordered B values", y="standardized introgressed ancestry"); plotD

# Plot recombination rate vs. introgressed ancestry for hominins and baboons
plotF <- ggplot() + 
  geom_smooth(data = hominin, aes(x=rank(recombination)/length(recombination), y=mean_ancestry_centered_sd), color="#366B7D", size=1.5, linetype = "longdash", method="lm", se=FALSE) + geom_smooth(data = baboon, aes(x=rank(rcr)/length(rcr), y=mean_ancestry_centered_sd), color="#366B7D", size=1.5,  method="lm", se=FALSE) + theme_classic() +
  theme(text=element_text(size=14), axis.text = element_text(color="black"), axis.text.x=element_blank(),axis.ticks.x=element_blank(), legend.position = "none", axis.title.x = element_text(vjust=-2), plot.margin = margin(l=20,b=10))  + 
  labs(x="rank-ordered recombination rate", y=""); plotF

# Plot all panels in figure 2
(plotA + theme(plot.margin = unit(c(0,0,20,0), "pt"))) + plotB + (plotC + theme(plot.margin = unit(c(0,0,20,0), "pt"))) + plotD + plotE + plotF + plot_layout(widths = c(1.3, 1))

ggsave("fig2.png")

#############################################################################################################################
# Figure 3
#############################################################################################################################

#############################################################################################################################
# Figure 3D-F
#############################################################################################################################

# Use the same baboon data frame from figure 2
# For right-hand panels, mean center anubis ancestry to 0 within recent and historic ancestry data to facilitate comparison of slopes
baboon$recent_ancestry_centered <- baboon$recent_ancestry - mean(baboon$recent_ancestry)
baboon$historical_ancestry_centered <- baboon$historical_ancestry - mean(baboon$historical_ancestry)

# We'll color recent ancestry blue/purple and historical ancestry brown 

# Plot fixed differences vs. introgressed recent and historical anubis ancestry
baboon$fixed[baboon$fixed > 30] <- 30 # classify all fixed differences greater than 30 in a "30+" bin which we'll note on the x-axis

plotDleft <- ggplot(data = baboon) +  geom_jitter(aes(x=fixed, y=recent_ancestry), alpha = 0.08, col="#5B507A") + geom_smooth(aes(x=fixed, y=recent_ancestry), method="lm", color="#5B507A", size=1.5, linetype = "longdash") + geom_jitter(aes(x=fixed, y=historical_ancestry), alpha = 0.08, col="#AD7748") + geom_smooth(aes(x=fixed, y=historical_ancestry), method="lm", color="#AD7748", size=1.5, linetype = "longdash") + theme_classic() + theme(text=element_text(size=25), axis.text = element_text(color="black")) + scale_x_continuous(labels=c("0","10","20","30+")) + xlab(label="") + scale_y_continuous(name= "proportion of anubis ancestry"); plotDleft

plotDright <- ggplot(data = baboon) + geom_smooth(aes(x=fixed, y=recent_ancestry_centered), method="lm", se=FALSE, color="#5B507A", size=1.5, linetype = "longdash") + geom_smooth(aes(x=fixed, y=historical_ancestry_centered), method="lm", se=FALSE, color="#AD7748", size=1.5, linetype = "longdash") + theme_classic() + theme(text=element_text(size=25), axis.text = element_text(color="black"))  +  coord_cartesian(ylim=c(-0.04, 0.04)) + scale_x_continuous(labels=c("0","10","20","30+")) + xlab(label="") + scale_y_continuous(name = "proportion of anubis ancestry\n- proportion of mean anubis ancestry"); plotDright 

# Plot ranked ordered B values vs. introgressed recent and historical anubis ancestry
plotEleft <- ggplot(data = baboon) +  geom_point(aes(x=rank(B), y=recent_ancestry), alpha = 0.08, col="#5B507A") +geom_smooth(aes(x=rank(B), y=recent_ancestry), method="lm", color="#5B507A", size=1.5, linetype = "longdash") + geom_point(aes(x=rank(B), y=historical_ancestry), alpha = 0.08, col="#AD7748") + geom_smooth(aes(x=rank(B), y=historical_ancestry), method="lm", color="#AD7748", size=1.5, linetype = "longdash") + theme_classic() + theme(text=element_text(size=25), axis.text = element_text(color="black"), axis.text.x=element_text(color="white"), axis.ticks.x=element_blank()) + xlab(label="") + ylab(label = "proportion of anubis ancestry"); plotEleft

plotEright <- ggplot(data = baboon) + geom_smooth(aes(x=rank(B), y=recent_ancestry_centered), method="lm", se=FALSE, color="#5B507A", size=1.5, linetype = "longdash") + geom_smooth(aes(x=rank(B), y=historical_ancestry_centered), method="lm", se=FALSE, color="#AD7748", size=1.5, linetype = "longdash") + theme_classic() + theme(text=element_text(size=25), axis.text = element_text(color="black"), axis.text.x=element_text(color="white"),axis.ticks.x=element_blank()) + xlab(label="") + ylab(label = "proportion of anubis ancestry\n- proportion of mean anubis ancestry") +  coord_cartesian(ylim=c(-0.04, 0.04)); plotEright

# Plot ranked ordered recombination rate vs. introgressed recent and historical anubis ancestry
plotFleft <- ggplot(data = baboon) +  geom_point(aes(x=rank(rcr), y=recent_ancestry), alpha = 0.08, col="#5B507A") + geom_smooth(aes(x=rank(rcr), y=recent_ancestry), method="lm", color="#5B507A", size=1.5, linetype = "longdash") + geom_point(aes(x=rank(rcr), y=historical_ancestry), alpha = 0.08, col="#AD7748") + geom_smooth(aes(x=rank(rcr), y=historical_ancestry), method="lm", color="#AD7748", size=1.5, linetype = "longdash") + theme_classic() + theme(text=element_text(size=25), axis.text = element_text(color="black"), axis.text.x=element_text(color="white"),axis.ticks.x=element_blank()) + xlab(label="") + ylab(label = "proportion of anubis ancestry"); plotFleft

plotFright <- ggplot(data = baboon)  + geom_smooth(aes(x=rank(rcr), y=recent_ancestry_centered), method="lm", se=FALSE, color="#5B507A", size=1.5, linetype = "longdash") + geom_smooth(aes(x=rank(rcr), y=historical_ancestry_centered), method="lm", se=FALSE, color="#AD7748", size=1.5, linetype = "longdash") + theme_classic() + theme(text=element_text(size=25), axis.text = element_text(color="black"), axis.text.x=element_text(color="white"), axis.ticks.x=element_blank()) + xlab(label="") + ylab(label = "proportion of anubis ancestry\n- proportion of mean anubis ancestry") + coord_cartesian(ylim=c(-0.04, 0.04)); plotFright

# Plot all panels in figures 3D-F
(((plotDleft + theme(plot.margin = unit(c(0,0,50,0), "pt"))) | plotDright) + xlab(label = "number of fixed differences per Mb") + theme(axis.title.x = element_text( hjust=2.25, vjust=-0.3))) / (((plotEleft + theme(plot.margin = unit(c(0,0,50,0), "pt"))) | plotEright) + xlab(label = "rank-ordered B values") + theme(axis.title.x = element_text(hjust=-6,vjust=-0.3))) / (((plotFleft + theme(plot.margin = unit(c(0,0,50,0), "pt"))) | plotFright) + xlab(label = "rank-ordered recombination rate") + theme(axis.title.x = element_text(hjust=3, vjust=-0.3)))

ggsave("fig3DEF.png") # this plot, I manually exported as an image in RStudio and save as "png" with 1000 width x 2000 height
