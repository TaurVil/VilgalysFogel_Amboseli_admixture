# Script for recreating figure 4A and 4B
# Note that the cartoon gene expression drawings beneath 4A and cartoon distribution in 4B were created in Keynote 8.1

# Load R libraries
library(ggplot2)

#############################################################################################################################
# Figure 4A
#############################################################################################################################

# Load Figure 4A results
gene <- read.table("fig4_example_gene.txt", header=T) # example gene MRPL2
# TAURAS FILE HERE # for q-q plot of p-values

# main panel
ggplot(gene) + 
  geom_violin(aes(x=as.factor(nearest_state), y=expression, fill=as.factor(nearest_state)), alpha=0.7, size=0.75) + 
  geom_boxplot(aes(x=as.factor(nearest_state), y=expression), width=0.1, size=0.75) + 
  scale_fill_manual(values=c("#FFED4F", "orange2", "#009E73")) + 
  scale_x_discrete(name="", labels=c("0" = "homozygous\nyellow", "1" = "heterozygous", "2" = "homozygous\nanubis")) + scale_y_continuous(name="normalized, residual\ngene expression") +
  theme_classic() + theme(legend.position = "none", axis.text = element_text(colour = "black"), text=element_text(size=20))

ggsave("fig4A_main.png")

# inset
# TAURAS CODE HERE

#############################################################################################################################
# Figure 4B
#############################################################################################################################

# Load Figure 4B results
results <- read.table("fig4BC_final_data/bootstrap_rho_results.txt", header=T)

# change 95% CI to standard error
results$high_se_DE <- results$rho_DE + results$sd_DE
results$low_se_DE <- results$rho_DE - results$sd_DE
results$high_se_nDE <- results$rho_nDE + results$sd_nDE
results$low_se_nDE <- results$rho_nDE - results$sd_nDE

ggplot(data=results) +
  geom_pointrange(aes(x=quantile, y=rho_DE, ymin=low_se_DE, ymax=high_se_DE), fill="#366B7D", color="#366B7D",alpha=0.8, size=1, shape=21, stroke=1.5) +
  geom_pointrange(aes(x=quantile+0.01, y=rho_nDE, ymin=low_se_nDE, ymax=high_se_nDE), fill="#5C438A", color="#5C438A",alpha=0.8, size=1, shape=21, stroke=1.5) +
  xlab("proportion of genes") + ylab("Spearman's rho") +
  coord_cartesian(ylim=c(0,0.2)) +
  theme_classic() +  theme(text=element_text(size=20), axis.text = element_text(color="black")) + xlab(label = "proportion of genes") # asterisks denoting significance were added later in keynote

# quantiles to add asterisks to in keynote
subset(results, bootstrap_pval<0.05)$quantile
# 0.15 0.20

ggsave("fig4B.png")
