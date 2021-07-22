# VilgalysFogel_Amboseli_admixture
## Selection against introgression and gene regulatory divergence in Amboseli baboons

This repository contains materials used to analyze the data from, and recreate the results in, Vilgalys, Fogel, et al. currently under review. These scripts include: 

### 01 Genotype calling: Scripts for trimming, mapping, and calling genotype data. 

### 02 Local ancestry calling: Local ancestry calling, including verifying our approach through simulation and quality control measures. 

### 03 Genomic resources for baboons: Annotation of the baboon genome, expanding upon Batra et al. 2020. 
In particular, this folder includes code and references for a liftover file between the human genome (hg19) and baboon genome (Panubis1.0), recombination rates (using data from Robinson et al. 2019), a measure of background selection (B values: McVicker et al. 2009), the annotation of CpG islands, and a liftover of human enhancers to the baboon genome. 

### 04 Structure of the baboon hybrid zone: Analyses and figures related to the extent of admixture beyond the Amboseli hybrid zone. 

### 05 Selection against introgression in Amboseli: Analyses and figured related to selection against anubis baboon ancestry in Amboseli baboons. 

### 06 Selection against regulatory divergence: Analyses and figures related to ancestry-associated expression, and selection against those regions.  

### 07 Predicting the genomic landscape of introgression: Analyses and figures related to predicting the extent of anubis introgression. 
There are two sub-sections here: (i) machine learning to identify genomic features predicting the extent of introgression and (ii) longitudinal analysis of changes in ancestry over time across the genome. 

Finally, we include here a RData file (VilgalysFogel_main_data_file.250kb.RData) which contains estimates of mean anubis ancestry for 250kb windows and genes across the baboon genome. This data file can be produced from the code and information available in sections 01-03, is produced by code in section 05, and is used for analyses in sections 05-07. 

For additional assistance, please contact Tauras Vilgalys (taur.vil@gmail.com) and Arielle Fogel (afogel29@gmail.com). 
