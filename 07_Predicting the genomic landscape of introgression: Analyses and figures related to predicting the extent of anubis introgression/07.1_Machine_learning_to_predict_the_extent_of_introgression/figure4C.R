# Script for recreating figure 4C

# Load R libraries
library(ggplot2)
library(tidyr)

#############################################################################################################################
# Figure 4C
#############################################################################################################################

# Load 4C results
results <- read.delim("glmnet_results.txt", header=T)

# Convert data from wide to long format for easier plotting
results2 <- results %>% pivot_longer(everything(), names_to=c("predictor"), values_to=c("estimate"))
# Should now have a data frame that is 4*200 rows (4 predictors, 200 estimates per predictor)
nrow(results2)==4*200 # TRUE

# Check to make sure nothing went wrong in the data frame conversion
summary(results$number.of.SNVs)==summary(subset(results2, predictor=="number.of.SNVs")$estimate) # all TRUE
summary(results$recombination.rate)==summary(subset(results2, predictor=="recombination.rate")$estimate) # all TRUE
summary(results$highly.differentiated.sites)==summary(subset(results2, predictor=="highly.differentiated.sites")$estimate) # all TRUE
summary(results$B)==summary(subset(results2, predictor=="B")$estimate) # all TRUE

ggplot(data = results2, aes(x=reorder(predictor, estimate), y=estimate)) + 
  geom_hline(yintercept = 0, color="grey", linetype="longdash", size=1.5) + 
  geom_jitter(fill="grey60", color="grey50", size=3, alpha=0.3, shape=21, stroke=1.5) +
  geom_boxplot(width=0.3, size=0.75, outlier.color = NA) +
  scale_x_discrete(labels=c("number.of.SNVs" = "number of SNVs", "recombination.rate" = "recombination rate", "highly.differentiated.sites" = "highly differentiated sites", "B" = "B statistic")) + scale_y_continuous(name="standardized\neffect size") +
  theme_classic() + theme(text=element_text(size=20), axis.text = element_text(color="black"), axis.title.y = element_blank()) +
  coord_flip()

ggsave("fig4C.png")
