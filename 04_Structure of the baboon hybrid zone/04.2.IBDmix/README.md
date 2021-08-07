## Estimate identity by descent (IBD) between anubis and yellow baboons. 

Results from LCLAE suggest some populations of yellow baboons previously thought to be unadmixed, specifically baboons used to found the SNPRC colony, may actual contain considerably anubis ancestry. However, these analyses rely upon reference panels of supposedly unadmixed individuals and admixture may extend further into the yellow baboon reference panel. Therefore, we used IBDmix (Chen et al. 2020, https://doi.org/10.1016/j.cell.2020.01.012) to estimate the amount of shared ancestry between yellow and anubis baboon individuals of different populations. IBDmix is a reference-free approach originally developed to detect archaic Neanderthal ancestry in modern African populations, meaning it does not assume any individuals are unadmixed. Specifically, IBDmix identifies sections of the genome which are IBD between individuals within a test population and a single individual representing the source of introgression. Because IBD is bidirectional (meaning it is unable to differentiate between anubis ancestry in yellow baboons and yellow ancestry in anubis baboons), we test all possible pairs of yellow and anubis baboons. We look for signatures of IBD that are shared across many potential source indivdiuals, reasoning that introgressed ancestry is unlikely to rise to high allele frequencies outside of the commonly recognized hybrid zone. 

A second limitation to IBD-based inference of introgression is that a baseline level of IBD is expected due to incomplete lineage sorting. To understand the expected rates of IBD in yellow an anubis baboons, we also apply IBDmix to estimate shared ancestry between yellow baboons with Guinea and hamadryas baboons, and between anubis baboons with yellow and Kinda baboons. These six species diverged approximately 1.5 million years ago into two clades: the northern clade of Guinea, anubis, and hamadryas baboons, and the southern clade of yellow, Kinda, and chacma baboons. As these other species are unlikely to have hybridized outside of their clade (due to strong geographic separation), estimates of IBD between different species serve as a proxy for the expected amount of IBD under incomplete lineage sorting alone without admixture. 


### Installation of IBDmix
```console 
# Local installation of up to date cmake 
./configure  --prefix=~/Programs/cmake/
gmake; gmake install 
# add cmake to the path
export PATH=$PATH:~/Programs/cmake/bin/

# install ibdmix: https://github.com/PrincetonUniversity/IBDmix
git clone https://github.com/PrincetonUniversity/IBDmix.git
cd IBDmix; mkdir build; cd build; cmake ..; cmake --build

# generates the `generate_gt` and `ibdmix` executables
~/Programs/IBDmix/IBDmix-master/build/src/
```

### Run IBDmix
We run IBDmix in parallel across chromosomes. For each run, we first grab genotype calls from previously called files (see Section 1) for yellow, anubis, and other source individuals. Genotype calls are then integrated for each data set using `generate_gt`, a function within IBDmix. IBDmix is then run sequentially for each source file and the results are exported as a series of text files.  

```console 
## have "source" files which contain the source individuals we want to iterate calling local ancestry over: 00_yellow_sources.list, 00_anubis_sources.list
## for each chromosome, run "run_IBDmix.sh"
sbatch --array=1-20 --mem=16G run_IBDmix.sh
```

### Estimate the mean proportion of the genome IBD between each test individual and source population
Following suggestions in Chen et al., we filter for tracts of IBD at least 50kb in length and with a LOD score greater than 10 in order to estimate the proportion of the genome shared between yellow and anubis baboons. For each baboon, we then estimate the proportion of the genome which is IBD with each possible source indivdual, and summarize these data for each baboon as the mean amount of IBD with each potential source population. These data are used in Fig 1C and Supplementary Tables, and are stored as `ibdmix_anubis_estimates.txt` and `ibdmix_yellow_estimates.txt`. 



```console
./01_mean_IBD.R
```

### Identify tracks of introgressed ancestry 
Next we sought to leverage IBDmix and LCLAE to remove putatively introgressed ancestry from our reference panels. To begin, we search for tracts of ancestry that are IBD with many putatively unadmixed source individuals (SNPRCanubis for yellow baboons; Mikumi for anubis baboons), reasoning that these sites are likely to be due to recent introgression events. 

Next we sought to identify the overlap between shared ancestry identified using LCLAE and IBDmix. 

.  When looking for overlap between multiple source individuals we applied more liberal thresholds for calling a region of the genome IBD, specifically tracts longer than 1000bp and LOD score greater than 4. While these criteria likely retain false positives, 


### For comparison to the SNPRC "yellow" founders, estimate IBD between Amboseli and all baboon species
We repeated the above process for nine high coverage Amboseli baboons to serve as a positive control of what , using all baboon species profiled as possible sources of ancestry. The files to run IBDmix are included here as `00_amboseli.list`, `00_amboseli_sources.list`, and `run_IBDmix_amboseli.sh`. Allele sharing between Amboseli baboons and each population was estimating using the follow R code. Unlike above, we don't collect where in the genome these tracts are located but rather simply estimate a mean proportion of the genome which is IBD between each Amboseli baboon and other baboon species. These data are summarized in `ibdmix_amboseli_estimates.txt` and used in Fig 1C. 

```console
./01b_mean_IBD_amboseli.R
```
