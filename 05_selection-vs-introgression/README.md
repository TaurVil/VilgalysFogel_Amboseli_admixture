## Selection against anubis introgression in the Amboseli baboons

To investigate whether signatures of selection against archaic introgression in humans are paralleled in baboons, we performed three tests. First, we asked whether regions of the genome that exhibit greater divergence between the parent species have reduced introgressed (anubis) ancestry in Amboseli. In humans, loci that contain fixed or near-fixed differences between humans and Neanderthals are more likely to be depleted of archaic ancestry (Vernot & Akey 2014 _Science_). Second, we tested whether regions of the genome that are predicted to be more affected by background selection also exhibit reduced introgressed ancestry, as reported for Neanderthal introgression into the human genome by Sankararaman et al. 2014 _Nature_. Third, we tested whether introgressed ancestry is positively correlated with recombination rate as seen in humans and other taxa (Schumer et al. 2018 _Science_).

For all three analyses, we calculated mean anubis ancestry for non-overlapping windows of the genome (100-1000 kb) in (i) all baboons from Amboseli (n=442); (ii) historical hybrids only (n=214); and (iii) recent hybrids only (n=188). In the main text, we report results using local ancestry estimates from LCLAE and 250 kb windows. However, the results are qualitatively unchanged if we use ancestry assignments based on Ancestry HMM (Table S3) or using alternative window sizes.

#### For each window size, get mean ancestry and genomic features

We will do this in two parts. First, we will calculate each feature for non-overlapping 25 kb windows, then we will combine these 25 kb windows into larger ones for later analyses. This is done to minimize the amount of computational time calculating each feature, while maintaining the flexibility to scale to larger window sizes. 

```console 
## get non-overlapping 25 kb windows across the genome and genomic features for them (we'll scale these smaller windows up to larger sizes for analyses)
## requires information on ancestry, recombination, B values, and the number of variable sites per window. 
## ancestry will be from amboseli_LCLAE_tracts.txt, available on Zenodo. With this code, we exclude excludes tracts less than 1 kb. 
## B values and recombination are called from files available in Section 03. Information on variable sites was too large to be included on GitHub, but may be generated from the code available in Section 03 and is provided for each window 25kb window in `windows.25kb.RData`. 

./run.01.get_25kb_windows.R
## produces `windows.25kb.RData`, which has for each window, the genomic features (features) which includes columns for the mean ancestry of all Amboseli baboons, of recently admixed baboons, and of baboons containing only historical admixture. Information about recent ancestry is included as a copy of table S1 (ids) and a matrix with the mean ancestry per recently admixed individual and chromosome (recent) is also included. 
## this code can also be modified to print out the ancestry per individual and window. This isn't included because of memory limitations, but can be done by retaining `ancestry_per_individual`. This can then be scaled to arbitrary window sizes using the code appended at the end of `run.02.scale_to_larger_windows.R`. 
```

Scale up from 25 kb windows to larger window sizes that are multiples of 25. 

```console
./run.02.scale_to_larger_windows.R 
## adjust the disance `distance` and output distance name `d_name` at the top of the script
```

#### Comparison with human-Neanderthal data

To provide context for our findings, we compared our results from baboons to evidence for selection against archaic introgression in humans, focusing specifically on admixture with Neanderthals. All three of the tests we used have been previously reported in the literature (Vernot & Akey 2014 _Science_; Sankararaman et al. 2014 _Nature_; Schumer et al. 2018 _Science_), but using different human reference genomes and differing amounts of information (based on what was available at the time) on Neanderthal-human divergence and introgression. We therefore re-ran these tests here, which allowed us to follow an analysis pipeline that was as parallel to our analysis of the baboons as possible.

```console
## create matrix with genomic and ancestry information for 250 kb windows along the human genome. Processing these data involves previously published data files (the sources of these data are described in detail in the Supplementary Methods). 
## these data are saved in `windows.human.250kb.RData`
./run.03.human_data.R

## get results for the human data which recapitulate those previously reported
./run.03b.human_results.R

```

#### Test for associations between genomic features and mean ancestry

For each window size, test for signatures of selection between introgressed ancestry and the number of fixed differences between yellow and anubis baboons, B values, and recombination rate. Although they can also be recapitulated from scripts `run.02.scale_to_larger_windows.R` and `run.04.results.R`, we also include a specific script `run.04b.results_250kb.R` which recapitulates results presented in the main text using `VilgalysFogel_main_data_file.250kb_windows.RData`. 

```console 
./run.04.results.R
./run.04b.results_250kb.R
## results for the main text, SI, and plot for supplementary figure S7 (note Fig S7 in the text uses 250 kb windows)

./figures2_3DEF.R
## plotting the results per 250 kb window for Fig 2 and Fig 3D-F
## this uses the main data file, which can also be exported by 02.scale_to_larger_windows.R with a window size of 250 kb, and the file of human data from `windows.human.250kb.RData`.
```
