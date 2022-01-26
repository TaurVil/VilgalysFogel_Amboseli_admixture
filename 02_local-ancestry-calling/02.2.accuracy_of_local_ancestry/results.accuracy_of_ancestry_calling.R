## SI: ACCURRACY OF ANCESTRY CALLING RESULTS
library(data.table); library(ggplot2); library(reshape2)

load("./DATA/ACCURACY_OF_ANCESTRY_CALLING.RData")

for (i in 1:25) {
  agreement$agreement[agreement$name == unique(agreement$name)[i] & 
                        agreement$cov == '10x' & agreement$method=="LCLAE"] <- 
    windowsize$agreement[windowsize$name == unique(agreement$name)[i] & 
                           windowsize$cov == '10x' & windowsize$distance == 35]
  agreement$agreement[agreement$name == unique(agreement$name)[i] & 
                        agreement$cov == '1x' & agreement$method=="LCLAE"] <- 
    windowsize$agreement[windowsize$name == unique(agreement$name)[i] & 
                           windowsize$cov == '1x' & windowsize$distance == 35]
}

## Plot comparing methods
A <- ggplot(agreement, aes(x=method, y=agreement, color=cov)) + 
  geom_boxplot(outlier.shape=NA, show.legend=F, aes(color=factor(cov))) + 
  scale_color_manual("coverage", values = c('darkorange4','darkorange')) +
  xlab("ancestry calling method") +
  geom_jitter(position=position_jitterdodge(jitter.width = 0.3), stroke=0,size=3.5, alpha=0.5, aes(color=factor(cov))) + 
  theme_classic() + 
  theme(legend.position = c(0.1,.9), legend.text = element_text(size=14), legend.title = element_text(size=14)) + 
  ylab("proportion of the genome called correctly") +
  guides(colour = guide_legend(override.aes = list(size=6))) + 
  theme(axis.text = element_text(color="black"), text=element_text(size=16))

## Plot comparing window sizes 
B <- ggplot(windowsize, aes(x=as.factor(distance), y=agreement, color=cov)) + 
  geom_boxplot(outlier.shape=NA, show.legend=F) + 
  scale_color_manual("coverage", values = c('darkorange4','darkorange')) +
  geom_jitter(position=position_jitterdodge(jitter.width = 0.2), stroke=0,size=2.5, alpha=0.5, aes(color=factor(cov))) + 
  theme_classic() + theme(legend.position = c(0.8,.2), legend.text = element_text(size=14), legend.title = element_text(size=14)) + 
  ylab("proportion of the genome called correctly") + xlab("window size (kb)") +
  guides(colour = guide_legend(override.aes = list(size=6))) + 
  theme(axis.text = element_text(color="black"), text=element_text(size=16))

## LCLAE: accuracy
  mean(agreement$agreement[agreement$cov == '10x' & agreement$method == "LCLAE"])
  sd(agreement$agreement[agreement$cov == '10x' & agreement$method == "LCLAE"])
  
  mean(agreement$agreement[agreement$cov == '1x' & agreement$method == "LCLAE"])
  sd(agreement$agreement[agreement$cov == '1x' & agreement$method == "LCLAE"])
  
  min(agreement$agreement[agreement$method == "LCLAE"])
  max(agreement$agreement[agreement$method == "LCLAE"])


## Accuracy with simulated reference panel
mean(simpanel$agreement)

## Unused plot comparing minimum delta
ggplot(minimumdelta, aes(x=as.factor(delta), y=agreement, color=cov)) + geom_boxplot(outlier.shape=NA) + 
  scale_color_manual("coverage", values = c('purple4', 'red')) +
  geom_jitter(position=position_jitterdodge(jitter.width = 0.1), size=2, aes(color=factor(cov))) + 
  theme_classic() + theme(legend.position = c(0.2,.15), legend.text = element_text(size=14), legend.title = element_text(size=14)) + 
  ylab("proportion of the genome called correctly") + xlab("minimum allele frequency difference")

## HMM: accuracy

  mean(agreement$agreement[agreement$cov == '10x' & agreement$method == "Ancestry-HMM"])
  sd(agreement$agreement[agreement$cov == '10x' & agreement$method == "Ancestry-HMM"])
  min(agreement$agreement[agreement$cov == '10x' & agreement$method == "Ancestry-HMM"])
  
  mean(agreement$agreement[agreement$cov == '1x' & agreement$method == "Ancestry-HMM"])
  sd(agreement$agreement[agreement$cov == '1x' & agreement$method == "Ancestry-HMM"])
  
  min(agreement$agreement[agreement$method == "Ancestry-HMM"])
  max(agreement$agreement[agreement$method == "Ancestry-HMM"])


## ADLIBS: accuracy

  mean(agreement$agreement[agreement$cov == '10x' & agreement$method == "ADLIBS"])
  sd(agreement$agreement[agreement$cov == '10x' & agreement$method == "ADLIBS"])
  min(agreement$agreement[agreement$cov == '10x' & agreement$method == "ADLIBS"])
  
  mean(agreement$agreement[agreement$cov == '1x' & agreement$method == "ADLIBS"])
  sd(agreement$agreement[agreement$cov == '1x' & agreement$method == "ADLIBS"])
  
  min(agreement$agreement[agreement$method == "ADLIBS"])
  max(agreement$agreement[agreement$method == "ADLIBS"])


### Greater accuracy in longer tracts 
subset(all_tracts, all_tracts$cov == "1x" & all_tracts$length > 50000) -> tmp
sum(tmp$correct)/sum(tmp$length)
subset(all_tracts, all_tracts$cov == "1x" & all_tracts$length < 1000) -> tmp
sum(tmp$correct)/sum(tmp$length)

subset(all_tracts, all_tracts$cov == "10x" & all_tracts$length > 50000) -> tmp
sum(tmp$correct)/sum(tmp$length)
subset(all_tracts, all_tracts$cov == "10x" & all_tracts$length < 1000) -> tmp
sum(tmp$correct)/sum(tmp$length)
rm(tmp)


## Plot agreement by tract length 
library(tidyquant); library(patchwork)
res <- subset(all_tracts, all_tracts$length > 0 & all_tracts$cov == "10x")
res <- res[order(res$length),]
lengths1 <- ggplot(res, aes(y=correct/length, x=log10(length))) + geom_point(alpha=0.08, stroke=0, col='dimgrey') + 
  geom_ma(n=70, ratio=0, color='darkorange4', linetype = "solid") + theme_classic() + 
  xlab("log10(called tract length)") + ylab("proportion of track called correctly") + ggtitle("high coverage (10x)") + 
  theme(axis.text = element_text(color="black"), text=element_text(size=16))

res <- subset(all_tracts, all_tracts$length > 0 & all_tracts$cov == "1x")
res <- res[order(res$length),]
lengths2 <- ggplot(res, aes(y=correct/length, x=log10(length))) + geom_point(alpha=0.08, stroke=0, col='dimgrey') + 
  geom_ma(n=70, ratio=0, color='darkorange', linetype = "solid") + theme_classic() + 
  xlab("log10(called tract length)") + ylab("proportion of track called correctly") + ggtitle("low coverage (1x)") + 
  theme(axis.text = element_text(color="black"), text=element_text(size=16))
#lengths1 + lengths2


A + (B / (lengths1 + lengths2))

mean(res$length); median(res$length)

