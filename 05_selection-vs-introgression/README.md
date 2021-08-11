## Selection against anubis introgression in Amboseli baboons

To investigate whether signatures of selection against archaic introgression in humans are paralleled in baboons, we performed three tests. First, we asked whether regions of the genome that exhibit greater divergence between the parent species have reduced introgressed (anubis) ancestry in Amboseli. In humans, loci that contain fixed or near-fixed differences between humans and Neanderthals are more likely to be depleted of archaic ancestry (Vernot & Akey 2014). Second, we tested whether regions of the genome that are predicted to be more affected by background selection also exhibit reduced introgressed ancestry, as reported for Neanderthal introgression into the human genome by Sankararaman et al. 2014. Third, we tested whether introgressed ancestry is depleted in regions of the genome with low mean recombination rate, as expected if deleterious alleles are selectively eliminated along with linked sequence and observed in humans (Schumer et al. 2018).

For all three analyses, we calculated mean anubis ancestry for non-overlapping windows of the genome (100-1000 kb) in (i) all baboons from Amboseli (n=442); (ii) historical hybrids only (n=214); and (iii) recent hybrids only (n=188). We report results from local ancestry calling with LCLAE and 250 kb windows in the main text. However, the results are qualitatively unchanged if we use ancestry assignments based on Ancestry HMM (Table S3) or using alternate window sizes.

#### For each window size, get mean ancestry and genomic features

```console 
## get non-overlapping 25 kb windows across the genome, and genomic features for them. We'll scale these smaller windows up to our larger sizes for analyses in each window size. 
## requires information for ancestry, recombination, B values, and the number of variable sites per window. 
## ancestry will be from amboseli_LCLAE_tracts.txt, availalbe on Zenodo, and excludes tracts less than 1 kb. 
## B values and recombination are called from files available in Section 03. Information of variable sites was too large to be included on GitHub, but may be generated from the code available in Section 03 and is provided for each window in the produced RData file. 

./r01.get_25kb_windows.R
## produces `windows.25kb.RData`, which has for each window the genomic features (features), the ancestry calls per individuals (ancestry_per_individual), and information about recent ancestry divided between a copy of table S1 (ids) and a matrix with the mean ancestry per recent individual and chromosome (recent). 







```