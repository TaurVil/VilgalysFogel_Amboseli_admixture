## Simulations and accuracy of local ancestry calling

Before calling local ancestry using LCLAE (Section 02.1), we sought to identify the best method to call local ancestry from low coverage baboon genotype data. To do so, we simulated local ancestry data for 25 individuals. We then simulated sequencing data for these individuals, based on their ancestry, and called local ancestry based on 10x or 1x sequencing coverage using three programs suitable for low coverage data. 


### Simulate local ancestry

We simulated local ancestry using SELAM (https://github.com/russcd/SELAM) to produce a population with a similar proportion of anubis ancestry and similar tract lengths to high coverage individuals in the Amboseli population. We then simulated high and low coverage sequencing data for these individuals and called genotypes following the procedure described in Section 01 with minor changes for a different computing cluster (exact code provided here). 

First, we will do some setup involving getting a reduced genome for the last 4 baboon chromosomes and getting the header for the vcf file. To do so, manually run through `run.00.setup.sh`. 

Next, we'll use SELAM to simulate ancestry along 4 autosomes (discarding an additional X chromosome). The output file must be ordered in increasing generations, and the output will apend to rather than overwrite an existing output file. 

```console
~/Programs/SELAM/src/SELAM -d ./DATA/demography.txt -o ./DATA/output.txt --seed 42 -c 5 77 84 63 63 1
# Demography: a population of 1000 individuals, 70% yellow with 7% and 3% new immigrants per generation
# Samples: we'll output 5 individuals at 5, 7, 10, 15, and 20 generations. Because calling local ancestry for is harder for shorter tracts and even 10 generations after admixture (in this simulation) produce shorter tracts than in high coverage Amboseli baboons, these simulations likely underestimate our ability to accurately call ancestry within Amboseli. 
#-c says to call 4 chromosomes, with the lengths given in morgans, which correspond to the smallest chromosomes in the baboon genome (17-20).
# The last chromosome is inherited from only the maternal line, hence why it has a weird output, and we exclude it from downstream analysis. 
deactivate
```

For each individual output by SELAM, let's remove the comment lines and the last chromosome which is only inherited along the maternal lineage. Then let's convert the output file from two haplotypes into ancestry calls of 0/1/2 anubis alleles along the chromosome. We'll call these the "true" tracts because they are known from the simulation. 

```console
ls i*output.txt | sed 's/.output.txt//g' > 00names 
for f in `ls i*output.txt`; do grep -v '^#' $f | grep -v -P "\t4\t" > tmp; mv tmp $f; done

mkdir true_tracts; module load R; for f in `cat 00names`; do sed -e s/NAME/$f/g run.01.SELAM_to_tracts.R > s.$f.sh; sbatch s.$f.sh; rm s.$f.sh; done
## There are checks in place to make sure each call is biallelic and length > 0
mv i*.txt true_tracts/; mv slurm-* true_tracts/ ## Cleanup
```

Having simulated ancestry tracts, we will now get genotypes for each individual. We'll do so by randomly drawing alleles based on an individual's ancestry state at the locus and the mean allele frequencies of the SNPRC anubis and yellow baoon colony founders. After generating genotypes for each individual, we will simulate 10x coverage sequencing using NEAT-genreads (Zachary Stephens, https://github.com/zstephens/neat-genreads and https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5125660/), map to the Panubis1 genome, and downsample to 1x and 2x coverage. 

```console
# get allele frequencies for yellow and anubis individuals that make up the reference panel
# apply to each chromosome we simulated data for (e.g. chr17-20)
# this is done using vcftools' `freq` command: vcftools --gzvcf masked_yellow_and_anubis.vcf.gz --keep 00_SNPRCanubis.list --chr CHROM --freq --out CHROM.anubis

# get genotype calls, by simulated individual and chromosome
mkdir simulated_vcfs; module load R; for f in `cat 00names`; do for g in `cat 00chroms`; do cat run.02.create_sample_vcf.R | sed -e s/NAME/$f/g | sed -e s/SCAF/$g/g > g.$f.$g.R; sbatch g.$f.$g.R; rm g.$f.$g.R; done; done 

# generate 10x coverage in SE, 100 bp reads for each individual and chromosome
mkdir sim_reads; for f in `cat 00names`; do for g in `cat 00chroms`; do sed s/INDIV/$f/g run.03.get_fasta.sh | sed s/CHROMO/$g/g > g.$f.sh; sbatch g.$f.sh; rm g.$f.sh; done; done

# combine chromosomes and map reads using bowtie2
mkdir mapped_bams; sbatch --array=1-25 --mem=8G run.04.map_sim_reads.sh

# downsample each sample to 2x and 1x coverage (so we have 10x, 2x, and 1x coverage data for each sample)
sbatch --array=1-25 --mem=6G run.05.downsample.sh
```

With bam files for each individual at each coverage, we will now call genotypes using GATK. We'll create one vcf file for each coverage and merge these with genotype calls for unadmixed yellow and anubis individuals (as was done for the full dataset including Amboseli individuals). 

```console
# call genotypes for each set of bams
cd mapped_bams/; ls i*bam | sed 's/.bam//g' > ../00_bams.txt; cd ..
mkdir gVCF; sbatch --array=1-75 --mem=16G run.06.gvcf.sh

# perform joint genotype calling for each set of bams 
sbatch --array=1-3 --mem=16G run.07.merge_gvcfs.sh
```

To control for sampling error in the reference panel, we also resampled a reference population of yellow and anubis baboons with allele frequencies drawn from the "true" frequencies observed from the actual reference panel. These individuals were simulated by creating "local ancestry tracts" that were entirely anubis or yellow, simulating a vcf for each sample using `run.02.create_sample_vcf.R`, simulating coverage, and calling genotypes based on the actual coverage of reference panel individuals. Using this alternate reference panel had little effect on ancestry calling and therefore is not included here.  

### Call local ancestry for simulated individuals

At this point, we have genotypes for each coverage titled `tmp.COVERAGE.vcf.gz`. We will now call local ancestry for these simulated individuals. 

To do so, we used LCLAE (Wall et al. 2016), ADLIBS (Schaefer et al. 2017), and Ancestry HMM (Corbett-Detig et al. 2017). Out of these three methods, LCLAE performed the best for low coverage data and is therefore the focus of this section. Our pipelines for running Ancestry HMM and ADLIBS are described separately in `running_adlibs.md` and `running_ancestryhmm.md`. Nevertheless, we summarize the results of these different analyses and, in the next section of this README, provide data and code necessary to recapitulate Fig. S2. 

**LCLAE (Wall et al. 2016)** uses genotype likelihoods extracted from a merged vcf file containing the reference and unadmixed individuals. 

```console
# get anubis and yellow genotypes for chromosomes 17-20 from the filtered files in Section 01. Combine into `refpanel.vcf.gz`
# check chromosome names for each file and make sure both either have a `chr` prefix before each chromosome number (e.g. chr1, chr2) or lack one (e.g. 1, 2)

# merge test samples with reference panel genotypes
# get genotype likelihoods for LCLAE
sbatch --array=1-3 --mem=16G run.08.get_genotypelikelihoods.sh

# use vcftools to get the depth (idepth), relatedness, and individual names from one of the CommonCalls vcfs. These help verify sample order for the *.h files and confirm sequencing depth of simulated data after mapping is similar to what was targetted 

# create anubis.h and yellow.h which are the number of individuals for each reference population (24 and 7, respectively) followed by their position in the vcf file as a space separated file (e.g. `24 1 2 14 15 ...`, if the anubis individuals were samples 1, 2, 14, and 15). 

# create directory to output results
mkdir raw_calls

# run LCLAE for each coverage (COVERAGE should be replaced with 1x, 2x, 10x, etc.)
# the array should specify the "test" individuals in the VCF
# this file calls `02_AIMS.txt` and `02_windows.txt` which cover a range of LCLAE parameters we test over. For the main results, we used values of 0.2 and 35 kb. 
for f in `cut -f 1 01_targetted_chroms.bed`; do sed -e s/CHROMOSOME/$f/g run.09.run_lclae.COVERAGE.sh > tmp.sh; sbatch --array=32-56 --mem=5G tmp.sh; done 
# array is run for samples 32-56 in the vcf file, which are our 25 test individuals simulated above but not the reference panel samples
# make sure that each file completes. LCLAE sometimes has random seg fault errors, but it works to rerun the script until all files have text. Alternately, files and conditions can be run individually which may work better for the last couple iterations to finish. 

# this produces one file per chromosome per individual, with information in the title regarding what parameters were used. 
# attach chromosome labels, then collapse chromosomes into single file 
for name in `cat test_samples.txt`; do for d in `cat 02_AIMs.txt`; do for window in `cat 02_windows.txt`; do for cov in `cat 01_sets.txt`; do for chrom in `cat 01_targetted_chroms.bed | cut -f 1`; do sed -i 's/^/\t/' ./raw_calls/sw.$window.$d.$name.$cov.$chrom.txt; sed -i s/^/$chrom/ ./raw_calls/sw.$window.$d.$name.$cov.$chrom.txt; done; cat ./raw_calls/sw.$window.$d.$name.$cov.chr* > ./raw_calls/sw.$window.$d.$name.$cov.txt; echo $name $d $window $cov; done; done; done; done 
rm raw_calls/*chr*.txt

# get LCLAE tracts after applying majority rule. 
# for majority rule, use the same distance as above except here it is distance from the focal SNP rather than window size. We tried modifying this separately but found it had little effect on the accuracy of ancestry calls. 

module load R; for name in `cat test_samples.txt`; do  for d in `cat 02_AIMs.txt`; do for window in `cat 02_windows.txt`; do for cov in `cat 01_sets.txt`; do sed -e s/INDIV/$name/g -e  s/DIFFERENCE/$d/g -e  s/DISTANCE/$window/g -e  s/COVERAGE/$cov/g run.07.calls_to_MajRule_tracts.R > r.$cov.$d.$window.R; sbatch --mem=20G --nice=50 r.$cov.$d.$window.R; sleep 10; rm r.$cov.$d.$window.R; fi; done; done; done; done 

# we now have ancestry tracts for each individual but have excluded the first and last 50 kb of each chromosome. 
# min50 and n20 in the title refer to a Majority Rule of at least 50% of sites and a minimum of 20 sites within the window to make a reliable call. Altering these parameters had little effect on the accuracy of ancestry calling. 
# tracts have not yet been filtered for length

```

### Accuracy of ancestry calls from simulated data

To evaluate how well different methods of ancestry calling performed, we compared the "true", simulated local ancestry tracts to those called using LCLAE, Ancestry HMM, and ADLIBS. After finding that LCLAE performed best for 1x coverage data, we then compared LCLAE results using different parameter choices to optimize ancestry calling for low coverage baboon data. 

The main component of this is an R function `agree` that compares `tract1` and `tract2`, and reports a matrix with the number of rows in tract1. For each row, `agree` returns for `tract1` the length of the tract, the ancestry state, and proportion of the tract with ancestry matching tract2. It also returns the number of base pairs assigned to each ancestry state in `tract2` for that tract. This function is saved in `./DATA/agree_function.RData`. The results for all methods and LCLAE parameters are summarized in `./DATA/ACCURACY_OF_ANCESTRY_CALLING.RData`. 

```console
## For each individual, read in the true ancestry tracts. 
#### make sure these tracts are biallelic

## then read in called local ancestry
## run agree and save the result to a growing file with information for all individuals and comparisons
res <- agree(test_tracts, true)
res$percent <- res$correct/res$length; res$name <- nom; res$cov <- tmp_cov; res$method <- tmp_method

## for LCLAE, retain the results for Fig S2C-D
## rbind(all_tracts, res) -> all_tracts

## report a summarized version for each individual and comparison
tmp_res <- as.data.frame(matrix(ncol=4, nrow=1)); colnames(tmp_res) <- c("name", "cov", "method", "agreement")
tmp_res$agreement <- sum(res$correct)/sum(res$length); tmp_res[,1] <- nom; tmp_res[,2] <- tmp_cov; tmp_res[,3] <- tmp_method
rbind(agreement, tmp_res) -> agreement
## Repeat for each individual and method, adding to `agreement` each time. 
```

Summarized data for the accuracy of local ancestry calling and code to recreate values in the text and Fig S2 can be found in `./results.accuracy_of_ancestry_calling.R` and `./DATA/ACCURACY_OF_ANCESTRY_CALLING.RData`. 
