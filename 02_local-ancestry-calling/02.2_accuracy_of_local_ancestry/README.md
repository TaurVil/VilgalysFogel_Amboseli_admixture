## Simulations and accuracy of local ancestry calling

Before calling genotypes using LCLAE (Section 02.1), we sought to identify the best method to call local ancestry from low coverage baboon genotype data. To do so, we simulated local ancestry data for 25 individuals. We then simulated sequencing data for these individuals and called local ancestry based on 10x or 1x sequencing coverage using three programs suitable for low coverage data. 


#### Simulate local ancestry

We simulated local ancestry using SELAM (https://github.com/russcd/SELAM) to produce a population with a similar proportion of anubis ancestry and similar tract lengths to high coverage individuals in the Amboseli population. We then simulated high and low coverage sequencing data for these individuals and called genotypes following the procedure described in Section 01 with minor changes for a different computing cluster (exact code provided here). 

First we will do some setup involving getting a reduced genome for the last 4 baboon chromosomes and getting the header for the vcf file. To do so, manually run through `run.00.setup.sh`. 

Next, we'll use SELAM to simulate ancestry along 4 autosomes (discarding an additional X chromosome). The output file must be ordered in increasing generations, and the output will apend to rather than overwrite an existing output file. 

```console
~/Programs/SELAM/src/SELAM -d ./DATA/demography.txt -o ./DATA/output.txt --seed 42 -c 5 77 84 63 63 1
# Demography: a population of 1000 individuals, 70% yellow with 7% and 3% new immigrants per generation
# Samples: we'll output 5 individuals at  5, 7, 10, 15, and 20 generations. Even 10 generations seems long for the mean Amboseli tract lengths (at least for the high coverage samples I compared them to, so this problem is actually harder than the actually problem posed in the real data)
#-c says to call 4 chromosomes, with the lengths given in morgans, which correspond to the smallest chromosomes in the baboon genome (17-20).
# The last chromosome is inherited from only the maternal line, hence why it has a weird output, and we exclude it from downstream analysis. 
deactivate
```

For each individual output by SELAM, let's remove the comment lines and the last chromosome which is only inheritted along the maternal lineage. Then let's convert the output file from two haplotypes into ancestry calls of 0/1/2 anubis alleles along the chromosome. We'll call these the "true" tracts because they are known from the simulation. 

```console
ls i*output.txt | sed 's/.output.txt//g' > 00names 
for f in `ls i*output.txt`; do grep -v '^#' $f | grep -v -P "\t4\t" > tmp; mv tmp $f; done

mkdir true_tracts; module load R; for f in `cat 00names `; do sed -e s/NAME/$f/g run.01.SELAM_to_tracts.R > s.$f.sh; sbatch s.$f.sh; rm s.$f.sh; done
## There are checks in place to make sure each call is biallelic and length > 0
mv i*.txt true_tracts/; mv slurm-* true_tracts/ ## Cleanup
```

Having simulated ancestry tracts, we now get genotypes for each individual. We'll do so by randomly drawing alleles based on an individuals ancestry state at the locus and the mean SNPRC founder anubis and yellow baboon allele frequencies (Section 3). After generating genotypes for each individual, we will simulate 10x coverage sequencing using NEAT-genreads (Zachary Stephens, https://github.com/zstephens/neat-genreads and https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5125660/), map to the Panubis1 genome, and downsample to 1x and 2x coverage. 

```console
# get genotype calls
mkdir simulated_vcfs; module load R; for f in `cat 00names`; do for g in `cat 00chroms`; do cat run.02.create_sample_vcf.R | sed -e s/NAME/$f/g | sed -e s/SCAF/$g/g > g.$f.$g.R; sbatch g.$f.$g.R; rm g.$f.$g.R; done; done 

# Generate 10x coverage in SE, 100bp reads for each individual and chromosome
mkdir sim_reads; for f in `cat 00names`; do for g in `cat 00chroms`; do sed s/INDIV/$f/g run.03.get_fasta.sh | sed s/CHROMO/$g/g > g.$f.sh; sbatch g.$f.sh; rm g.$f.sh; done; done

## Combine chromosomes, and map reads using bowtie2
mkdir mapped_bams; sbatch --array=1-25 --mem=8G run.04.map_sim_reads.sh

## Downsample each sample to 2x and 1x (as well as the 10x we generated)
sbatch --array=1-25 --mem=6G run.05.downsample.sh
```

With bam files for each individual at each coverage, we now call genotypes using GATK. We'll create one vcf file for each coverage and merge these with genotype calls for unadmixed yellow and anubis individuals (as was done for the full dataset including Amboseli individuals). 

```console
## Call genotypes for each set of bams
cd mapped_bams/; ls i*bam | sed 's/.bam//g' > ../00_bams.txt; cd ..
mkdir gVCF; sbatch --array=1-75 --mem=16G  run.06.gvcf.sh

## Joint genotype calling for each set of bams 
sbatch --array=1-3 --mem=16G run.07.merge_gvcfs.sh
```


#### Call local ancestry for simulated individuals


**LCLAE (Wall et al. 2016)**
W

```console
## Merge test samples with reference panel genotypes
module load java/1.8.0_45-fasrc01
java -jar ~/Programs/GenomeAnalysisTK.jar -T CombineVariants -R $path_genome --variant ./filt.$coverage.recode.vcf.gz --variant ../unadmixed_individuals/refpanel.vcf.gz -o merged.$coverage.vcf.gz -genotypeMergeOptions UNIQUIFY -L 01_targetted_chroms.bed

java -jar ~/Programs/GenomeAnalysisTK.jar -T SelectVariants -R $path_genome -V merged.$coverage.vcf.gz -select 'set == "Intersection"' -o CommonCalls.$coverage.vcf

## Exctract genotype likelihoods


```
