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

### Merge genotypes for the high coverage samples we will be using
```console 
## have "source" files which contain the source individuals we want to iterate calling local ancestry over: 00_yellow_sources.list, 00_anubis_sources.list

```

### Merge genotypes for the high coverage samples we will be using
```console

```

### For comparison to the SNPRC "yellow" founders, estimate IBD between Amboseli and all baboon species
```console
## "source" file saved as 00_amboseli_sources.txt
## 00_amboseli.list contains the 9 high coverage Amboseli baboon genomes

## Get Amboseli genotypes from: ~/panubis1_genotypes/calls_merged/04.merged_shared.$chrom.vcf.gz
## Get source genotypes from: ~/panubis1_genotypes/high_coverage/02.unadmixed_highcov.$chrom.recode.vcf.gz

mkdir ./chrom_vcfs/ ; mkdir IBDmix_by_chrom ## chrom vcfs contains intermediaries while IBDmix_by_chrom contains the real output
sbatch --array=1-20 --mem=16G run_IBDmix.sh
```
