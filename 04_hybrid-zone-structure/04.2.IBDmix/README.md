## Estimate identity by descent (IBD) between anubis and yellow baboons. 

Results from PCA and LCLAE suggest some populations of yellow baboons previously thought to be unadmixed, specifically baboons used to found the SNPRC colony, may actually contain considerable anubis ancestry. However, these analyses rely upon reference panels of supposedly unadmixed individuals and admixture may extend further into the yellow baboon reference panel. Therefore, we used IBDmix (Chen et al. 2020 _Cell_, https://doi.org/10.1016/j.cell.2020.01.012) to estimate the amount of shared ancestry between yellow and anubis baboon individuals from different populations. IBDmix is a reference-free approach, meaning it does not assume _any_ individuals are unadmixed, and was originally developed to detect archaic Neanderthal ancestry in modern African populations. Specifically, IBDmix identifies sections of the genome which are IBD between individuals within a "test" population and a single individual representing the "source" of introgression. Because IBD is bidirectional (meaning it is unable to differentiate between anubis ancestry in yellow baboons and yellow ancestry in anubis baboons), we tested all possible pairs of yellow and anubis baboons. We looked for signatures of IBD that are shared across many potential source indivdiuals, reasoning that introgressed ancestry is unlikely to rise to high allele frequencies outside of the commonly recognized hybrid zone. 

A second limitation to IBD-based inference of introgression is that a baseline level of IBD is expected due to incomplete lineage sorting. To understand the expected rates of IBD in yellow and anubis baboons, we also apply IBDmix to estimate shared ancestry between yellow baboons with Guinea and hamadryas baboons, and between anubis baboons with yellow and Kinda baboons. These six species diverged approximately 1.5 million years ago into two clades: the northern clade of Guinea, anubis, and hamadryas baboons, and the southern clade of yellow, Kinda, and chacma baboons. As these other species are unlikely to have hybridized outside of their clade (due to large geographic separation), estimates of IBD between different species serve as a proxy for the expected amount of IBD under incomplete lineage sorting alone without admixture. 


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
mkdir ./IBDmix_by_chrom; mkdir ./chrom_vcfs
sbatch --array=1-20 --mem=16G run_IBDmix.sh
```

### Estimate the mean proportion of the genome IBD between each test individual and source population
Following suggestions in Chen et al., we filter for tracts of IBD at least 50 kb in length and with a LOD score greater than 10 in order to estimate the proportion of the genome shared between yellow and anubis baboons. For each baboon, we then estimate the proportion of the genome which is IBD with each possible source indivdual, and summarize these data for each baboon as the mean amount of IBD with each potential source population. These results are used in Fig 1C and Supplementary Tables, and are stored as `ibdmix_anubis_estimates.txt` and `ibdmix_yellow_estimates.txt`. 

```console
./r01_mean_IBD.R
```

### Identify tracks of introgressed ancestry 
Next, we sought to leverage both IBDmix and LCLAE to remove putatively introgressed ancestry from our reference panels. To begin, we search for tracts of ancestry that are IBD with many putatively unadmixed source individuals (SNPRCanubis for yellow baboons; Mikumi for anubis baboons), reasoning that these sites are likely to be due to recent gene flow. While identifying these regions, we used more liberal thresholds for calling regions of the genome IBD, specifically tracts larger than 1000 bp with a LOD score greater than 4. While these criteria likely retain some false positives, we reasoned that patterns of sharing across individuals would add additional confidence and that it would be more harmful to retain true cases of introgression than mask excess variation. 

We return tracts for each individual that are likely to contain heterospecific ancestry based on 30, 50, or 70% of possible source individuals. Zipped versions of these files are included in the `RESULTS` folder.

```console
## get the tracts that are IBD for each sample
./r02a_concensus_IBD_tracts.R
# saved as an RData object (`./IBDmix_tracts.RData`) and six text files (`./RESULTS/ibdmix_tracts.SPECIES.THRESHOLD.txt`), with species as either "anubis" or "yellow" and thresholds as either 30, 50, or 70% of possible source individuals. 
```

If IBD is a signature of gene flow, we expect that results from IBDmix should overlap with LCLAE. For LCLAE, we used local ancestry calls for unadmixed yellow and anubis individuals as described in Section 02.1. These files are included here in the `DATA/lclae` folder. These results reveal substantial overlap between IBDmix and LCLAE for SNPRC yellow baboons, far exceeding what is expected by chance. This script also includes code to produce Fig S5. 

```console
## integrate IBDmix results with LCLAE results to get the overlap between the two approaches, both genome-wide and at base pair resolution
./r02b_IBD_and_LCLAE.R
```

Based on the above results, we conclude that there is likely previously undetected gene flow between anubis and yellow baboons which may affect our ability to call local ancestry. We therefore mask regions of the genome which contain possible gene flow between our reference panel individuals, using the intersection of local ancestry results using LCLAE and IBDmix requiring 50% of source individuals. With these criteria we increase the density of ancestry informative markers, decrease the occurrence of pedigree inconsistencies, and rarely drop additional informative markers from our dataset. A masked vcf file (`masked_yellow_and_anubis.vcf.gz`) is available on Zenodo and produced in Section 01.3 from `yes_intersect_50.bed`, which we create here. We also confirmed results were highly consistent using masked and unmasked genotypes (>90% similarity in local ancestry, all major results agree), or alternate masking criteria such as not filtering LCLAE tracts or using more/less stringent IBDmix criteria. 

```console
## get file of tracts to mask

## for IBDmix, use ibdmix_tracts.anubis.50.txt and ibdmix_tracts.yellow.50.txt, which were produced by ./r02b_IBD_and_LCLAE.R and are included in the RESULTS directory
## for LCLAE, use "./RESULTS/lclae.txt" which can be produced using ./r02b_IBD_and_LCLAE.R and the LCLAE output files contained in the DATA directory

## attach 'chr' to the IBDmix tracts (which just have the chromosome number). LCLAE already has the `chr`
for f in `ls ./ibdmix_tracts.*`; do sed -i s/^/chr/g $f; sed -i s/chrchr/chr/g $f; echo $f; done 

## get subsets of LCLAE tracts for high coverage samples. The ./DATA/00_anu.list and ./DATA/00_yel.list refer to 28 and 18 individuals, respectively, which include high coverage data from SNPRC colony founders (n=7 yellow; 24 anubis), the Baboon Genome Diversity Project (anubis from Aberdare and SNPRC; 1 yellow from Mikumi), and Mikumi (n=10). 
## get subsets for low coverage samples too. This is the individuals not in the previous set, and includes 4 Mikumi baboons, 7 Maasai Mara baboons, and 6 WNPRC baboons. 
for f in `cat ./DATA/00_anu.list`; do grep $f masking_tracts/lclae.txt >> masking_tracts/lclae_highcov.txt; echo $f; done 
for f in `cat ./DATA/00_yel.list`; do grep $f masking_tracts/lclae.txt >> masking_tracts/lclae_highcov.txt; echo $f; done 
cp masking_tracts/lclae.txt masking_tracts/lclae_lowcov.txt; for f in `cat ./00_anu.list`; do grep -v $f masking_tracts/lclae_lowcov.txt > tmp.txt; mv tmp.txt masking_tracts/lclae_lowcov.txt; echo $f; done; for f in `cat ./00_yel.list`; do grep -v $f masking_tracts/lclae_lowcov.txt > tmp.txt; mv tmp.txt masking_tracts/lclae_lowcov.txt; echo $f; done 

## get intersect of IBDmix and LCLAE results for high coverage samples
cat masking_tracts/ibdmix_tracts.*.50.txt | grep $name > tmp.ibd.bed; grep $name masking_tracts/lclae_highcov.txt | sed -e 's/\.5//g' > tmp.lclae.bed; bedtools intersect -a tmp.ibd.bed -b tmp.lclae.bed | sed "s/$/\t${name}/"  >> masking_tracts/none_intersect_50.bed ; rm tmp.ibd.bed; rm tmp.lclae.bed 

## merge intersection with LCLAE for low coverage samples 
cat masking_tracts/none_intersect_50.bed masking_tracts/lclae_lowcov.txt | sed -e 's/\.5//g' > ./yes_intersect_50.bed
rm -r masking_tracts

```

### For comparison to the SNPRC "yellow" baboon founders, estimate IBD between Amboseli and all baboon species
We repeated the above analyses for nine high coverage Amboseli baboons to serve as a positive control of what , using all baboon species as possible sources of ancestry. The files to run IBDmix are included here as `00_amboseli.list`, `00_amboseli_sources.list`, and `run_IBDmix_amboseli.sh`. Allele sharing between Amboseli baboons and each population was estimating using the follow R code. Unlike above, we do not focus on where in the genome these tracts are located but rather simply estimate a mean proportion of the genome which is IBD between each Amboseli baboon and other baboon species. These data are summarized in `ibdmix_amboseli_estimates.txt` and used in Fig 1C. 

```console
./run_IBDmix_amboseli.sh

./01b_mean_IBD_amboseli.R
```

### For comparison to the SNPRC "yellow" baboon founders, estimate IBD between Amboseli and all baboon species
Figure 1C can be reproduced from the processed results in this section, based on the raw code and data which have been made available. 

```console
./figure1C.R
```
