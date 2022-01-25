# Script for recreating figure S7 PSMC plot

# Load R libraries
library(ggplot2)
library(dplyr)

# Plot PSMC results (years before present vs. estimated Ne)

load("PSMC_results.RData.", verbose=T)
#Loading objects:
#  Ne_estimates
#  time_estimates

Ne_estimates <- as.data.frame(Ne_estimates) 
time_estimates <- as.data.frame(time_estimates) 
Ne_estimates$row <- rownames(Ne_estimates)
time_estimates$row <- rownames(time_estimates)

Ne_estimates2 <-  pivot_longer(Ne_estimates, V1:V253)
names(Ne_estimates2)[3] <- "Ne"
time_estimates2 <-  pivot_longer(time_estimates, V1:V253)
names(time_estimates2)[3] <- "time"

all_data <- merge(Ne_estimates2, time_estimates2, by=c("row", "name"))

# Plot up to 2 million years ago and upt to 100000 Ne (to actually see all of the replicate samples)
ggplot(data=all_data) + geom_line(aes(time, Ne, group=name), color="black", alpha=0.1) + scale_x_continuous(name="years before present (millions)") + scale_y_continuous(name=expression("estimated N"[e])) + theme_classic() + theme(text=element_text(size=25), legend.position = "none", , axis.text = element_text(color="black"))+coord_cartesian(xlim=c(0,2), ylim=c(0,100000))
