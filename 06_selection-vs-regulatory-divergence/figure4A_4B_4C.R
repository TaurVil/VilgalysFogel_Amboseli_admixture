# Script for recreating figure 4A, 4B, and 4C
# Note that the cartoon gene expression drawings beneath 4A and cartoon distribution in 4C were created in Keynote 8.1

# Load R libraries
library(ggplot2)

#############################################################################################################################
# Figure 4A
#############################################################################################################################

# Load Figure 4A results

# main panel
gene <- read.table("./RESULTS/fig4_example_gene.txt", header=T) # example gene MRPL2

ggplot(gene) + 
  geom_violin(aes(x=as.factor(nearest_state), y=expression, fill=as.factor(nearest_state)), alpha=0.7, size=0.75) + 
  geom_boxplot(aes(x=as.factor(nearest_state), y=expression), width=0.1, size=0.75) + 
  scale_fill_manual(values=c("#FFED4F", "orange2", "#009E73")) + 
  scale_x_discrete(name="", labels=c("0" = "homozygous\nyellow", "1" = "heterozygous", "2" = "homozygous\nanubis")) + scale_y_continuous(name="normalized, residual\ngene expression") +
  theme_classic() + theme(legend.position = "none", axis.text = element_text(colour = "black"), text=element_text(size=20)) -> a

ggsave("fig4A_main.png")

# inset is produced based on a qq plot using all permutations. Code and data are not provided here due to memory constraints, but it is included when running r05.mash_linear_models.R

#############################################################################################################################
# Figure 4B
#############################################################################################################################
reduction <- read.table("./RESULTS/reduced_ancestry_in_DE_genes_per_individual.txt", header=T)
colnames(reduction) <- c('qqq', 't_data')

ggplot(data=reduction, aes(x=as.factor(qqq), y=t_data, group=qqq)) +
  geom_hline(yintercept=0, linetype="dashed", color="grey60") +
  geom_violin(alpha=0.7, size=0.75, fill="grey60") + 
  geom_boxplot(width=0.1, size=0.75) +
  scale_x_discrete(name="proportion of genes", 
                   labels=c("0.1" = "top 10%\nvs.\nbottom 10%", "0.15" = "15%", "0.2" = "top 20%\nvs.\nbottom 20%", "0.25" = "25%","0.3" = "top 30%\nvs.\nbottom 30%", "0.35" = "35%","0.4" = "top 40%\nvs.\nbottom 40%")) + 
  scale_y_continuous(name="\u0394 anubis ancestry\n(least DE - most DE)", limits=c(-0.13, 0.13))  +  
  theme_classic() +  theme(text=element_text(size=20), axis.text = element_text(color="black"), axis.text.x = element_text(face="italic", size=12), axis.title.x=element_text(vjust=-1))  -> b 

ggsave("fig4B.png")

#############################
# Figure 4C
#############################################################################################################################

# Load Figure 4C results
results <- read.table("./RESULTS/bootstrap_rho_results.txt", header=T)

# Quantiles to add asterisks to
subset(results, bootstrap_pval<0.05)$quantile
# 0.15 0.20

# Plot rho plus/minus SD
ggplot(data=results) +
  geom_pointrange(aes(x=quantile, y=rho_DE, ymin=rho_DE - sd_DE, ymax=rho_DE + sd_DE), fill="#0192CB", color="#0192CB", size=1, shape=21, stroke=1.5) +
  geom_pointrange(aes(x=quantile+0.01, y=rho_nDE, ymin=rho_nDE -  sd_nDE, ymax=rho_nDE + sd_nDE), fill="#965F8E", color="#965F8E", size=1, shape=21, stroke=1.5) +
  scale_x_continuous(name="proportion of genes", labels=c("0.1" = "top 10%\nvs.\nbottom 10%", "0.15" = "15%", "0.2" = "top 20%\nvs.\nbottom 20%", "0.25" = "25%","0.3" = "top 30%\nvs.\nbottom 30%", "0.35" = "35%","0.4" = "top 40%\nvs.\nbottom 40%"), breaks=seq(0.1,0.4, 0.05)) + scale_y_continuous("Spearman's rho", limits=c(0, 0.2)) +  theme_classic() +  theme(text=element_text(size=20), axis.text = element_text(color="black"), axis.text.x = element_text(face="italic", size=12), axis.title.x=element_text(vjust=-1)) + annotate("text", x = 0.15, y = 0.18, label = "*", size=8) +  annotate("text", x = 0.2, y = 0.18, label = "*", size=8) -> c

ggsave("fig4C.png")
