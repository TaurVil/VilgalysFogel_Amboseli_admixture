## Overlap between IBDmix and LCLAE
library(IRanges); library(data.table); library(ggplot2); library(patchwork)

load("./RESULTS/IBDmix_tracts.RData")

## Read in LCLAE names
lclae_names <- read.delim("./DATA/01_lclae_samples.list", header=F)

#### Pull in lclae results for each sample ############
lclae <- NULL; for (n in lclae_names$V1) {
  nom <- paste("./DATA/full_refpanel.lclae35kb.p2.majrule35kb.min50.n20.exclude50kb.",n,".tracts.txt", sep="")
  tmp <- fread(nom)
  tmp$name <- n
  rbind(lclae, tmp) -> lclae
}; rm(nom, tmp, n)

yel_sw <- unlist(read.delim("./DATA/00_yel_SW.list", header=F))
yel_mikumi <- unlist(read.delim("./DATA/00_yel_wall.list", header=F))
anu_sw <- unlist(read.delim("./DATA/00_anu_SW.list", header=F))
anu_mara <- unlist(read.delim("./DATA/00_anu_mara.list", header=F))
anu_tulane <- unlist(read.delim("./DATA/00_anu_tulane.list", header=F))

lclae$pop <- "bgp_anubis"; lclae$pop[lclae$name %in% yel_sw] <- "SW yellow"; lclae$pop[lclae$name %in% yel_mikumi] <- "Mik yellow"
lclae$pop[lclae$name %in% anu_sw] <- "SW anubis"; lclae$pop[lclae$name %in% anu_mara] <- "Mara anubis"; lclae$pop[lclae$name %in% anu_tulane] <- "Tulane anubis"

lclae$spec <- "anubis"; lclae$spec[lclae$name %in% yel_sw] <- "yellow"; lclae$spec[lclae$name %in% yel_mikumi] <- "yellow"
lclae$spec[lclae$name %in% anu_sw] <- "anubis"; lclae$spec[lclae$name %in% anu_mara] <- "anubis"; lclae$spec[lclae$name %in% anu_tulane] <- "anubis"

#### print out lclae tracts that don't match the species ############
lclae_wrong_tracts <- lclae
lclae_wrong_tracts <- subset(lclae_wrong_tracts, (lclae_wrong_tracts$spec == 'anubis' & lclae_wrong_tracts$state < 2) | (lclae_wrong_tracts$spec == 'yellow' & lclae_wrong_tracts$state > 0))
write.table(lclae_wrong_tracts, "./RESULTS/lclae.txt", row.names=F, col.names=T, sep="\t", quote=F)

#### calculate mean ancestry per individual based on LCLAE ############
lclae_res <- as.data.frame(matrix(ncol=3, nrow=nrow(lclae_names))); colnames(lclae_res) <- c("name", "prop_anubis", "prop_wrong")
for (i in 1:nrow(lclae_res)) {
  lclae_res[i,1] <- as.character(lclae_names$V1[i])
  tmp <- subset(lclae, lclae$name == lclae_names$V1[i] & lclae$length > 1000)
  tmp$length <- tmp$end - tmp$start
  lclae_res[i,2] <- sum(tmp$length * tmp$state)/(2*sum(tmp$length))
}; rm(i, tmp)

lclae_res$pop <- "1_Aberdares anubis"; lclae_res$pop[lclae_res$name %in% yel_sw] <- "6_Southwest yellow"; lclae_res$pop[lclae_res$name %in% yel_mikumi] <- "5_Mikumi yellow"
lclae_res$pop[lclae_res$name %in% c(anu_sw, "panu_L142", "panu_LIV5")] <- "3_Southwest anubis"; lclae_res$pop[lclae_res$name %in% anu_mara] <- "2_Masai Mara anubis"; lclae_res$pop[lclae_res$name %in% anu_tulane] <- "4_Tulane anubis"

lclae_res$spec <- "anubis"; lclae_res$spec[lclae_res$name %in% yel_sw] <- "yellow"; lclae_res$spec[lclae_res$name %in% yel_mikumi] <- "yellow"
lclae_res$spec[lclae_res$name %in% anu_sw] <- "anubis"; lclae_res$spec[lclae_res$name %in% anu_mara] <- "anubis"; lclae_res$spec[lclae_res$name %in% anu_tulane] <- "anubis"

lclae_res$cov <- 'hi'
lclae_res$cov[! (lclae_res$name %in% yel$name | lclae_res$name %in% anu$name)] <- 'low'
lclae_res$colchoice <- 'purple4'; lclae_res$colchoice[lclae_res$cov == 'low'] <- 'red'

factor(lclae_res$pop) -> lclae_res$x

dev.off()

set.seed(1); ggplot(data=lclae_res, aes(x=x, y=prop_anubis)) + 
  scale_color_manual("coverage", values = c('#366B7D','#9BC6D4'), labels=c("high (>10x)", "low (<5x)")) +
  geom_boxplot(outlier.shape=NA, show.legend=F, size=1) + 
  geom_jitter(width = 0.4, aes(col=cov), size=4, alpha=0.75, stroke=0) + 
  theme_classic() + xlab("") + ylab("estimated proportion anubis ancestry") + 
  theme(legend.position = c(0.1,.15), legend.text = element_text(size=14), legend.title = element_text(size=14)) + 
  guides(colour = guide_legend(override.aes = list(size=6))) + 
  theme(axis.text = element_text(color="black"), text=element_text(size=16))
## Pdf version saved as "LCLAE_refpanels.pdf", 12x6, Fig S5A

## statistics for the text
mean(lclae_res$prop_anubis[lclae_res$pop %like% "Mikumi yellow" & lclae_res$cov == "low"])
sd(lclae_res$prop_anubis[lclae_res$pop %like% "Mikumi yellow" & lclae_res$cov == "low"])
mean(lclae_res$prop_anubis[lclae_res$pop %like% "Mikumi yellow" & lclae_res$cov == "hi"])
sd(lclae_res$prop_anubis[lclae_res$pop %like% "Mikumi yellow" & lclae_res$cov == "hi"])
mean(lclae_res$prop_anubis[lclae_res$pop %like% "Southwest yellow" & lclae_res$cov == "hi"])
sd(lclae_res$prop_anubis[lclae_res$pop %like% "Southwest yellow" & lclae_res$cov == "hi"])
1-lclae_res$prop_anubis[lclae_res$name %like% "panu_30877"]

for (i in unique(lclae_res$pop)) {
  print(i)
  print(mean(lclae_res$prop_anubis[lclae_res$pop == i]))
  print(median(lclae_res$prop_anubis[lclae_res$pop == i]))
  
}
print(mean(lclae_res$prop_anubis[lclae_res$pop %like% "Mikumi yellow" & lclae_res$cov == 'low']))
print(mean(lclae_res$prop_anubis[lclae_res$pop %like% "Mikumi yellow" & lclae_res$cov == 'hi']))


#### convert LCLAE to 1 for incorrect ancestry, 0 for correct ancestry, to match IBDmix ############

lclae$state2 <- lclae$state
lclae$state2[lclae$state == 2 & lclae$spec == 'yellow'] <- 1
lclae$state2[lclae$state == 1 & lclae$spec == 'anubis'] <- 0
lclae$state2[lclae$state == 2 & lclae$spec == 'anubis'] <- 1

#### get LCLAE information compared to IBDmix ####

#### Genome-wide comparison: Fig S5B-C ####
## length of states that are called incorrectly; tend to be much shorter than tracts assigned to the correct species
hist(log10(lclae$length[lclae$state == 1]), breaks=100); mean(lclae$length[lclae$state == 1], na.rm=T)
## 
for (n in yel$name) {
  subset(lclae, lclae$name == n & lclae$state == 1) -> tmp
  yel$lclae[yel$name == n] <- sum(tmp$length, na.rm=T)/sum(chroms$length)
}; rm(n, tmp)

plot(yel$lclae ~ rowMeans(yel[,2:29]), xlab="mean IBDmix", ylab="LCLAE", frame=F, col='blue', pch=16, ylim=c(0,.4)); abline(a=0,b=1,col='red')
summary(lm(yel$lclae ~ rowMeans(yel[,2:29]))); abline(lm(yel$lclae ~ rowMeans(yel[,2:29])))

## for each individual, get the proportion yellow/anubis IBD from IBDmix30,50,and 70
anu$ibd70 <- anu$ibd50 <- anu$ibd30 <- NA
for (i in 1:nrow(anu)) {
  tmp_anu <- anu$name[i]
  tmp_tracts <- subset(ibd_mix_tracts_30anu, ibd_mix_tracts_30anu$name == tmp_anu)
  anu$ibd30[i] <- sum(tmp_tracts$number * tmp_tracts$length)/sum(tmp_tracts$length)
  tmp_tracts <- subset(ibd_mix_tracts_50anu, ibd_mix_tracts_50anu$name == tmp_anu)
  anu$ibd50[i] <- sum(tmp_tracts$number * tmp_tracts$length)/sum(tmp_tracts$length)
  tmp_tracts <- subset(ibd_mix_tracts_70anu, ibd_mix_tracts_70anu$name == tmp_anu)
  anu$ibd70[i] <- sum(tmp_tracts$number * tmp_tracts$length)/sum(tmp_tracts$length)
}
yel$ibd70 <- yel$ibd50 <- yel$ibd30 <- NA
for (i in 1:nrow(yel)) {
  tmp_yel <- yel$name[i]
  tmp_tracts <- subset(ibd_mix_tracts_30yel, ibd_mix_tracts_30yel$name == tmp_yel)
  yel$ibd30[i] <- sum(tmp_tracts$number * tmp_tracts$length)/sum(tmp_tracts$length)
  tmp_tracts <- subset(ibd_mix_tracts_50yel, ibd_mix_tracts_50yel$name == tmp_yel)
  yel$ibd50[i] <- sum(tmp_tracts$number * tmp_tracts$length)/sum(tmp_tracts$length)
  tmp_tracts <- subset(ibd_mix_tracts_70yel, ibd_mix_tracts_70yel$name == tmp_yel)
  yel$ibd70[i] <- sum(tmp_tracts$number * tmp_tracts$length)/sum(tmp_tracts$length)
}
rm(i, tmp_anu, tmp_yel, tmp_tracts)
## add in proportion from LCLAE
for (n in yel$name) {
  subset(lclae, lclae$name == n & lclae$state == 1) -> tmp
  yel$lclae[yel$name == n] <- sum(tmp$length, na.rm=T)/sum(chroms$length)
}; rm(n, tmp)
for (n in anu$name) {
  subset(lclae, lclae$name == n & lclae$state == 1) -> tmp
  anu$lclae[anu$name == n] <- sum(tmp$length, na.rm=T)/sum(chroms$length)
}; rm(n, tmp)

to_plot <- as.data.frame(cbind(rep(yel$lclae,3), c(yel$ibd70, yel$ibd50, yel$ibd30), c(rep("ibd30",nrow(yel)), rep("ibd50",nrow(yel)), rep("ibd70",nrow(yel)))))
colnames(to_plot) <- c("lclae", "ibdmix", "threshold")
to_plot$lclae <- as.numeric(to_plot$lclae); to_plot$ibdmix <- as.numeric(to_plot$ibdmix)

yel_plot <- ggplot(to_plot, aes(x=ibdmix, y=lclae, col=threshold)) + 
  scale_color_manual("proportion of IBDmix\nsource individuals", values = c('skyblue4','skyblue2', 'skyblue'), labels=c("70%", "50%", "30%")) +
  geom_point(size=4, alpha=0.75, stroke=0) + 
  geom_abline(slope=1, intercept = 0, col='red', linetype='dashed') +
  theme_classic() + xlim(c(0,0.5)) + ylim(c(0,0.5)) +
  xlab("shared ancestry estimated using IBDmix") + 
  ylab("introgression estimated using LCLAE") +
  theme(legend.position = c(0.3,.9), legend.text = element_text(size=14), legend.title = element_text(size=14)) + 
  guides(colour = guide_legend(override.aes = list(size=6))) + 
  theme(legend.title.align=0.5, axis.text = element_text(color="black"), text=element_text(size=16))

to_plot <- as.data.frame(cbind(rep(anu$lclae,3), c(anu$ibd70, anu$ibd50, anu$ibd30), c(rep("ibd30",nrow(anu)), rep("ibd50",nrow(anu)), rep("ibd70",nrow(anu)))))
colnames(to_plot) <- c("lclae", "ibdmix", "threshold")
to_plot$lclae <- as.numeric(to_plot$lclae); to_plot$ibdmix <- as.numeric(to_plot$ibdmix)

anu_plot <- ggplot(to_plot, aes(x=ibdmix, y=lclae, col=threshold)) + 
  scale_color_manual("proportion of IBDmix\nsource individuals", values = c('skyblue4','skyblue2', 'skyblue'), labels=c("70%", "50%", "30%")) +
  geom_point(size=4, alpha=0.75, stroke=0, show.legend=F) + 
  geom_abline(slope=1, intercept = 0, col='red', linetype='dashed') +
  theme_classic() + xlim(c(0,0.35)) + ylim(c(0,0.35)) +
  xlab("shared ancestry estimated using IBDmix") + 
  ylab("introgression estimated using LCLAE") +
  theme(legend.position = c(0.8,.2), legend.text = element_text(size=14), legend.title = element_text(size=14)) + 
  guides(colour = guide_legend(override.aes = list(size=6))) + 
  theme(legend.title.align=0.5, axis.text = element_text(color="black"), text=element_text(size=16))
rm(to_plot)
yel_plot + anu_plot ## 12x6, ibdmix_vs_lclae.pdf  # Fig S5B and Fig S5C


#### site-specific comparison ####
## agree function defined in Section 02.2, used to compare LCLAE to each IBDmix threshold
agree <- function(tract1, tract2) {
  chroms <- unique(tract1$chrom)
  tract2$number -> tract2$state
  as.numeric(tract1$start) -> tract1$start
  as.numeric(tract1$end) -> tract1$end
  as.numeric(tract2$start) -> tract2$start
  as.numeric(tract2$end) -> tract2$end
  for (x in 1:length(chroms)) {
    #print("starting chrom")
    subset(tract1, tract1$chrom == chroms[x]) -> t1
    subset(tract2, tract2$chrom == chroms[x]) -> t2
    res <- as.data.frame(matrix(ncol=2, nrow=nrow(t1)))
    colnames(res) <- c('length', 'correct')
    res$chrom <- chroms[x]
    for (y in 1:nrow(t1)) {
      subset(t2, (t2$end > t1$start[y] & t2$start < t1$end[y])) -> t3
      t3$start[t3$start < t1$start[y]] <- t1$start[y]
      t3$end[t3$end > t1$end[y]] <- t1$end[y]
      t3$length <- as.numeric(t3$end) - as.numeric(t3$start)
      sum(t3$length) -> res$length[y]
      subset(t3, t3$state == t1$state[y]) -> t4
      sum(t4$length) -> res$correct[y]
    }
    if (x ==1) {res -> lengths} else {rbind(lengths,res) -> lengths}
    #print(paste("Done with",x,sep=" "))
  }
  return(lengths)
}

## create res_yel and res_anu, with columns referring to the proportion of LCLAE tracts called at a IBDmix threshold (e.g. P_lclae30) and the proportion of the genome called using IBDmix and LCLAE
## this code can be modified to call IBDmix as well, but this is very time consuming because IBD returns lots of small tracts
res <- as.data.frame(matrix(ncol=2, nrow=18)); colnames(res) <- c("name", "P_lclae30")
for (n in 1:length(yel$name)) {
  res$name[n] <- as.character(yel$name[n])
  tmp_l <- subset(lclae, lclae$name == yel$name[n] & lclae$state2 == 1) 
  tmp_l$state <- 1
  tmp_i <- subset(ibd_mix_tracts_30yel, ibd_mix_tracts_30yel$name == yel$name[n] & ibd_mix_tracts_30yel$number >= 1); tmp_i$chrom <- paste("chr",tmp_i$chrom,sep="")
  r_tmp <- agree(tmp_l, tmp_i)
  sum(r_tmp$correct[!is.na(r_tmp$length)])/sum(tmp_l$length[!is.na(tmp_l$length)]) -> res$P_lclae30[n]; 
  #r_tmp <- agree(tmp_i, tmp_l)
  #sum(r_tmp$correct[!is.na(r_tmp$length)])/sum(tmp_i$length[!is.na(tmp_i$length)]) -> res$P_ibdmix30[n]; print("done")
  res$ibd30[n] <- sum(tmp_i$length, na.rm=T)/sum(chroms$length)
  print(n)
  tmp_i <- subset(ibd_mix_tracts_50yel, ibd_mix_tracts_50yel$name == yel$name[n] & ibd_mix_tracts_50yel$number >= 1); tmp_i$chrom <- paste("chr",tmp_i$chrom,sep="")
  r_tmp <- agree(tmp_l, tmp_i)
  sum(r_tmp$correct[!is.na(r_tmp$length)])/sum(tmp_l$length[!is.na(tmp_l$length)]) -> res$P_lclae50[n]; 
  #r_tmp <- agree(tmp_i, tmp_l)
  #sum(r_tmp$correct[!is.na(r_tmp$length)])/sum(tmp_i$length[!is.na(tmp_i$length)]) -> res$P_ibdmix50[n]; print("done")
  res$ibd50[n] <- sum(tmp_i$length, na.rm=T)/sum(chroms$length)
  print(n)
  tmp_i <- subset(ibd_mix_tracts_70yel, ibd_mix_tracts_70yel$name == yel$name[n] & ibd_mix_tracts_70yel$number >= 1); tmp_i$chrom <- paste("chr",tmp_i$chrom,sep="")
  r_tmp <- agree(tmp_l, tmp_i)
  sum(r_tmp$correct[!is.na(r_tmp$length)])/sum(tmp_l$length[!is.na(tmp_l$length)]) -> res$P_lclae70[n]; 
  #r_tmp <- agree(tmp_i, tmp_l)
  #sum(r_tmp$correct[!is.na(r_tmp$length)])/sum(tmp_i$length[!is.na(tmp_i$length)]) -> res$P_ibdmix70[n]; print("done")
  res$ibd70[n] <- sum(tmp_i$length, na.rm=T)/sum(chroms$length)
  print(n)
  
  res$lclae[n] <- sum(tmp_l$length, na.rm=T)/sum(chroms$length)
  
}; res -> res_yel; rm(res, r_tmp, tmp_l, n, tmp_i)
res <- as.data.frame(matrix(ncol=2, nrow=28)); colnames(res) <- c("name", "P_lclae30")
for (n in 1:length(anu$name)) {
  res$name[n] <- as.character(anu$name[n])
  tmp_l <- subset(lclae, lclae$name == anu$name[n] & lclae$state2 == 0) 
  tmp_l$state <- 1
  tmp_i <- subset(ibd_mix_tracts_30anu, ibd_mix_tracts_30anu$name == anu$name[n] & ibd_mix_tracts_30anu$number >= 1); tmp_i$chrom <- paste("chr",tmp_i$chrom,sep="")
  r_tmp <- agree(tmp_l, tmp_i)
  sum(r_tmp$correct[!is.na(r_tmp$length)])/sum(tmp_l$length[!is.na(tmp_l$length)]) -> res$P_lclae30[n]; 
  #r_tmp <- agree(tmp_i, tmp_l)
  #sum(r_tmp$match[!is.na(r_tmp$length)])/sum(tmp_i$length[!is.na(tmp_i$length)]) -> res$P_ibdmix30[n]; print("done")
  res$ibd30[n] <- sum(tmp_i$length, na.rm=T)/sum(chroms$length)
  print(n)
  tmp_i <- subset(ibd_mix_tracts_50anu, ibd_mix_tracts_50anu$name == anu$name[n] & ibd_mix_tracts_50anu$number >= 1); tmp_i$chrom <- paste("chr",tmp_i$chrom,sep="")
  r_tmp <- agree(tmp_l, tmp_i)
  sum(r_tmp$correct[!is.na(r_tmp$length)])/sum(tmp_l$length[!is.na(tmp_l$length)]) -> res$P_lclae50[n]; 
  #r_tmp <- agree(tmp_i, tmp_l)
  #sum(r_tmp$match[!is.na(r_tmp$length)])/sum(tmp_i$length[!is.na(tmp_i$length)]) -> res$P_ibdmix50[n]; print("done")
  res$ibd50[n] <- sum(tmp_i$length, na.rm=T)/sum(chroms$length)
  print(n)
  tmp_i <- subset(ibd_mix_tracts_70anu, ibd_mix_tracts_70anu$name == anu$name[n] & ibd_mix_tracts_70anu$number >= 1); tmp_i$chrom <- paste("chr",tmp_i$chrom,sep="")
  r_tmp <- agree(tmp_l, tmp_i)
  sum(r_tmp$correct[!is.na(r_tmp$length)])/sum(tmp_l$length[!is.na(tmp_l$length)]) -> res$P_lclae70[n]; 
  #r_tmp <- agree(tmp_i, tmp_l)
  #sum(r_tmp$match[!is.na(r_tmp$length)])/sum(tmp_i$length[!is.na(tmp_i$length)]) -> res$P_ibdmix70[n]; print("done")
  res$ibd70[n] <- sum(tmp_i$length, na.rm=T)/sum(chroms$length)
  print(n)
  
  res$lclae[n] <- sum(tmp_l$length, na.rm=T)/sum(chroms$length)
  
}; res -> res_anu; rm(res, r_tmp, tmp_l, n, tmp_i)

## Over 75% of LCLAE calls are represented in the IBDmix 30 set
mean(res_yel$P_lclae30); mean(res_anu$P_lclae30)
## This drops to 50-60% for IBDmix with 70% of source individuals, but that is obviously a more stringent threshold
mean(res_yel$P_lclae70); mean(res_anu$P_lclae70)

## fisher exact test for enrichment of overlap between LCLAE and IBDmix 70 for anubis (the lowest overlap)
fisher.test(10000000*matrix(ncol=2, nrow=2, c(.5*0.0689, .5*0.0689, 0.076-(.5*0.0689), 1-.5*0.0689-0.076)))


