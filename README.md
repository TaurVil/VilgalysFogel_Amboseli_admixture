# VilgalysFogel_Amboseli_admixture
## Selection against admixture and gene regulatory divergence in a long-term primate field study

![alt text](https://github.com/TaurVil/VilgalysFogel_Amboseli_admixture/blob/main/baboon.jpg?raw=true)
<sup>*Young baboons playing in the Amboseli basin of southern Kenya*</sup>

This repository contains materials used to analyze the data from, and recreate the results in, Vilgalys, Fogel, et al. *bioRxiv*. These scripts are organizing in the following seven directories: 

1. **Genotype calling**: Trimming, mapping, and calling genotypes from whole genome resequencing data.

2. **Local ancestry calling**: Local ancestry estimation, including verifying our approach through simulations, quality control measures, and orthogonal approaches. 

3. **Genomic resources for baboons**: Annotation of the baboon genome, expanding upon data available from Batra et al. 2020 *GigaScience* (https://doi.org/10.1093/gigascience/giaa134). 

4. **Structure of the baboon hybrid zone**: Analyses and figures related to the geographical and temporal extent of admixture in and around the Amboseli hybrid population.

5. **Selection against introgression in Amboseli**: Analyses and figures related to selection against introgressed anubis ancestry in the Amboseli baboons. 

6. **Selection against regulatory divergence**: Analyses and figures related to ancestry-associated gene expression and selection against those regions of the genome.  

7. **Predicting the genomic landscape of introgression**: Analyses and figures related to predicting anubis introgression across the genome and over time. 

Finally, we include here a RData file (`VilgalysFogel_main_data_file.250kb.RData`) which contains estimates of mean anubis ancestry for 250 kb windows across the baboon genome. This file contains 4 tables and 5 other variables. The 4 tables are: (i) ids: a copy of Table S1 containing information for Amboseli baboons; (ii) recent: a matrix of individuals (rows) and chromsoomes (columns) containing the mean anubis ancestry per chromosome of individuals with a recent anubis ancestor known from the Amboseli pedigree; and 2 copies of genomic information and mean ancestry per 250kb window, the version of which titled `to_analyze` has been filtered for recombination rates higher than `max_rcr`. These files contain columns for the mean recombination rate (rcr), mean B values (B), and the number of fixed differences (fixed). This data file is produced in Section 05. 

Thank you for your interest in our work!

For additional assistance, please contact Tauras Vilgalys (taur.vil@gmail.com) and Arielle Fogel (afogel29@gmail.com). 
