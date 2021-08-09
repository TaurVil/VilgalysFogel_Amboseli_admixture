## Simulations and accuracy of local ancestry calling

Before calling genotypes using LCLAE (Section 02.1), we sought to identify the best method to call local ancestry from low coverage baboon genotype data. To do so, we simulated 


#### Simulating local ancestry

We simulated local ancestry using SELAM (https://github.com/russcd/SELAM) to produce a population with a similar proportion of anubis ancestry and similar tract lengths to high coverage individuals in the Amboseli population. We then randomly sampled anubis or yellow baboon genotypes (based on the simulated local ancestry and the allele frequency in each species) at each variable site to provide genotype for these individuals, before simulating 1x, 2x, and 10x coverage sequencing data. These data were mapped to the Panubis1 genome and used to call genotypes following the procedure described in Section 01 with minor changes for a different computing cluster (code provided here). 

```console
## First we will do some setup involving getting a reduced genome for the last 4 baboon chromosomes, getting the header for the vcf file

###Run SELAM
# must have output file in increasing generations
# must not have existing output files (will append to rather than overwrite)

#Load in modules and pre-requisites
module load gcc; module load samtools; module load python/2.7.6-fasrc01; module load virtualenv; module load gsl; 
# First time we will need to create the virtual environment: virtualenv venv/
source ~/local_ancestry/simulated_data/venv/bin/activate; pip install --upgrade pip==6.0.7
# and install SELAM: wget https://github.com/russcd/SELAM/archive/master.zip; unzip master.zip; cd ~/Programs/SELAM/src; make

#So let's output 4 chromosomes for each individual. 
~/Programs/SELAM/src/SELAM -d try1_demography.txt -o try1_output.txt --seed 42 -c 5 77 84 63 63 1
# Demography: a population of 1000 individuals, 70% yellow with 7% and 3% new immigrants per generation
# Samples: we'll output 5 individuals at  5, 7, 10, 15, and 20 generations. Even 10 generations seems long for the mean Amboseli tract lengths (at least for the high coverage samples I compared them to, so this problem is actually harder than the actually problem posed in the real data)
#-c says to call 4 chromosomes, with the lengths given in morgans, which correspond to the smallest chromosomes in the baboon genome (17-20).
# The last chromosome is inherited from only the maternal line, hence why it has a weird output, and we exclude it from downstream analysis. 
deactivate

ls i*output.txt | sed 's/.output.txt//g' > 00names 

## Let's clean up the tracks a little bit, removing the comment lines and the last chromosome (#4), which is from the maternal lineage only
for f in `ls i*output.txt`; do grep -v '^#' $f | grep -v -P "\t4\t" > tmp; mv tmp $f; done

##We must convert ancestry output from two haplotypes per chromosome into a 0/1/2 ancestry calls (the equivalent of the LCLAE output)
mkdir tracts; module load R; for f in `cat 00names `; do sed -e s/NAME/$f/g run.01.SELAM_to_tracts.R > s.$f.sh; sbatch s.$f.sh; rm s.$f.sh; done
## There are checks in place to make sure each call is biallelic and length > 0
mv i*.txt tracts/; mv slurm-* tracts/ ## Cleanup

# get a vcf file for each simulated sample, containing genotypes randomly drawn from the two populations based on ancestry state 
mkdir simulated_vcfs; module load R; for f in `cat 00names`; do for g in `cat 00chroms`; do cat run.02.create_sample_vcf.R | sed -e s/NAME/$f/g | sed -e s/SCAF/$g/g > g.$f.$g.R; sbatch g.$f.$g.R; rm g.$f.$g.R; done; done 

##Simulate reads using NEAT-genreads, a tool developed by Zachary Stephens, https://github.com/zstephens/neat-genreads and https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5125660/
# Generate 10x coverage in SE, 100bp reads
mkdir sim_reads; for f in `cat 00names`; do for g in `cat 00chroms`; do sed s/INDIV/$f/g run.03.get_fasta.sh | sed s/CHROMO/$g/g > g.$f.sh; sbatch g.$f.sh; rm g.$f.sh; done; done
mv tracts/ true_tracts/

## Combine chromosomes, and map reads using bowtie2
mkdir mapped_bams; sbatch --array=1-25 --mem=8G run.04.map_sim_reads.sh

## Downsample each sample to 2x and 1x (as well as the 10x we generated)
sbatch --array=1-25 --mem=6G run.05.downsample.sh

## Call genotypes for each set of bams
cd mapped_bams/; ls i*bam | sed 's/.bam//g' > ../00_bams.txt; cd ..
mkdir gVCF; sbatch --array=1-75 --mem=16G  run.06.gvcf.sh

## Joint genotype calling for each set of bams 
sbatch --array=1-3 --mem=16G run.07.merge_gvcfs.sh
## Clean out all pre-requisite files up to this point. 

```
