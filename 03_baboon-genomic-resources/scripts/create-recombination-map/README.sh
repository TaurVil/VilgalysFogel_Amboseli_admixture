## Get recombination maps for the panubis1 genome
# We'll be making an anubis map using the anubis SNPRC founder baboons (n=24)

## get genotype data from the raw, unfiltered vcf file 
## keep the 24 SNPRC high coverage anubis individuals, filter for sites with no missing data, that pass GATK hard filters, and have a minor allele frequency less than 5%. Also remove clusters of 3+ SNPs within 10 bp
## matches filtering criteria in Pfeifer 2020 (vervet monkey genetic map )

sbatch --array=1-20 --mem=12G ./run.01.prep_genotypes.sh 

## merge all SNPRC anubis genotypes into one file 
module load bcftools; bcftools concat ./SNPRCanubis.*.hwe.vcf.gz -O z -o ./SNPRCanubis.vcf.gz

## Phase with beagle
module load java; \java -Xmx4g -jar ~/Programs/beagle.08Jun17.d8b.jar gt=./SNPRCanubis.vcf.gz out=./SNPRCanubis_beagle ne=40000; tabix SNPRCanubis_beagle.vcf.gz

# Break up into individual chromosomes
for f in `seq 1 20`; do vcftools --gzvcf SNPRCanubis_beagle.vcf.gz --chr $f --recode --out $f.SNPRCanubis ; bgzip $f.SNPRCanubis.recode.vcf; tabix $f.SNPRCanubis.recode.vcf.gz; done
# Get individual-vcfs for each haplotype
for f in `seq 1 20`; do for g in `cat ./n24.anubis.SNPRCfounder.list`; do sed -e s/CHROM/$f/g run.01a.panubis1.sh | sed -e s/INDIV/$g/g > g.sh; sbatch g.sh; rm g.sh; done; done
# Run LDhelmet
for f in `seq 1 20`; do sed -e s/CHROM/$f/g run.02a.panubis1.sh > r.all.sh; sbatch --mem=65000 --nice r.all.sh; rm r.all.sh; done

