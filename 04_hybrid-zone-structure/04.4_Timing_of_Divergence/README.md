

## Estimate the timing of divergence between yellow and anubis baboons far from the hybrid zone

Applied to a F1 hybrid, pairwise sequential Markovian coalescence (PSMC: Li & Durbin, 2011 _Nature_) can be used to infer when gene flow ceased between two populations as this corresonds to a rapid increase in the effect population size to values approaching an infinitely large population (Cahill et al. 2016 _Proc B_). We therefore created pseudo-hybrid individuals by combining haploid genomes from a SNPRC anubis baboon and Mikumi yellow baboon, populations which to the best of our knowledge have not experienced admixture, and estimated the effective population size over time to confirm that the ancestors of yellow and anubis baboons ceased exchanging gene flow for some period of time. As expected, we observed a near-asymptotic increase in estimated effective population size, consistent with historical cessation of gene flow between species. 

```console 
# following README available at https://github.com/jacahill/hPSMC

# have bam files in this folder under "bams"

# install pu2fa from https://github.com/Paleogenomics/Chrom-Compare
# get a haploid genome file for each individual
mkdir haploid_fa
for f in `cat 00_yellow.list`; do sed -e s/MY_SAMPLE_NAME/$f/g run.get_haploid_files.sh > g.$f.sh; sbatch --mem=8G g.$f.sh; rm g.$f.sh; done 
for f in `cat 00_anubis.list`; do sed -e s/MY_SAMPLE_NAME/$f/g run.get_haploid_files.sh > g.$f.sh; sbatch --mem=8G g.$f.sh; rm g.$f.sh; done 

# get merged psmc fasta files 
mkdir pseudo_fasta
module load python/2.7.6-fasrc01
for f in `cat 00_yellow.list`; do for g in `cat 00_anubis.list`; do sed -e s/MY_YELLOW_NAME/$f/g -e s/MY_ANUBIS_NAME/$g/g run.get_pseudo_fasta.sh > g.$f.$g.sh; sbatch --mem=8G g.$f.$g.sh; rm g.$f.$g.sh; done ; done 

# download psmc: https://github.com/lh3/psmc
# run psmc on each pair of samples
for f in `cat 00_yellow.list`; do for g in `cat 00_anubis.list`; do sed -e s/MY_YELLOW_NAME/$f/g -e s/MY_ANUBIS_NAME/$g/g run.get_psmc.sh > g.$f.$g.sh; sbatch --mem=24G g.$f.$g.sh; rm g.$f.$g.sh; done ; done

# get output statistics scaled for generation time (11 years) and baboon mutation rate (0.57e-9 mutations per year: Wu et al. 2020, _PloS Biol_)
for f in `cat 00_yellow.list`; do for g in `cat 00_anubis.list | grep -v '1_'`; do python ../../Programs/hPSMC/PSMC_emit_last_iteration_coord.py -s10 -g11 -m0.0000000057 $f.$g.psmc > to_plot.p57.t11.$f.$g.txt; done; done 

```

Once all output files have been obtained, we will process them in R to get a single data object (see code below), included here as PSMC_results.RData. This data object can then be used to create figure S7 (`figureS7.R`).

```console
# in R

library(ggplot2)
# creat two matrixes, one for times and one for population sizes
Ne_estimates <- time_estimates <- NULL
for (i in list.files(pattern = "to_plot.p57.t11.")) {
  read.delim(i, sep=' ', header=F) -> d
  cbind(time_estimates, d[,1]) -> time_estimates
  cbind(Ne_estimates, d[,2]) -> Ne_estimates
  rm(d)
}; rm(i)


```


