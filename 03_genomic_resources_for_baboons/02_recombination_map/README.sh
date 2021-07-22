## Get recombination maps for the panubis1 genome
# We'll be making an anubis map using the SW founder baboons (n=24 SW anubis)

cd /data/tunglab/tpv/panubis1_genotypes/recombination

## Starting with the high coverage samples (in /data/tunglab/tpv/panubis1_genotypes/high_coverage and on github Amboseli_LocalAncestry/Resources/High Coverage Samples/)
module load bcftools; bcftools concat /data/tunglab/tpv/panubis1_genotypes/high_coverage/02.unadmixed_highcov.*.vcf.gz -O z -o ./unadmixed_highcov.vcf.gz

## Recombination map-specific genotype filtering for each dataset
# Filtering instructions from Pfeifer 2020 (vervet monkey genetic map)

## anubisSW
# Excluded sites with missing data to avoid errors and biases resulting from computational imputation
module load vcftools; vcftools --gzvcf unadmixed_highcov.vcf.gz --keep ../00_anu_SW.list --max-missing 1 --max-alleles 2 --recode --recode-INFO-all --out anubisSW --maf 0.05
## Remove clusters of SNPs within a 10bp window (GATK clusterSize=3, clusterWindowSize=10) and singletons within the target population
module load tabix; bgzip anubisSW.recode.vcf; module load java/1.8.0_45-fasrc01; module load tabix; tabix ./anubisSW.recode.vcf.gz
java -jar /data/tunglab/tpv/Programs/GenomeAnalysisTK.jar -T VariantFiltration -R ../Panubis1.nochromname.fa -V ./anubisSW.recode.vcf.gz -cluster 3 -window 10 -o ./anubisSW_cluster.vcf.gz
# Stringent filtering for false variants
## SNPs showing an excess of heterozygosity were removed. Specifically, a P-value for Hardy–Weinberg Equilibrium was calculated using the “—hardy” option in VCFtools v.0.1.13 (Danecek et al. 2011), and SNPs with P < 0.01 removed.
module load vcftools; vcftools --hardy --gzvcf anubisSW_cluster.vcf.gz --out anubisSW
module load R; R
library(data.table); fread("anubisSW.hwe") -> hwe; subset(hwe, hwe$P_HWE < 0.01) -> hwe; dim(hwe); write.table(hwe[,1:2], "anubisSW_hwe_to_remove.sites", row.names=F, col.names=F, sep="\t", quote=F)
vcftools --gzvcf anubisSW_cluster.vcf.gz --recode-INFO-all --remove-filtered-all --exclude-positions anubisSW_hwe_to_remove.sites --recode --out anubisSW_hwe
bgzip anubisSW_hwe.recode.vcf ; tabix anubisSW_hwe.recode.vcf.gz
## SNPs that could be recipricolly lifted over with the human genome? We'll ignore this, assuming the quality of the baboon genome is sufficient. 
## Remove fixed alleles: done earlier with the maf frequency

## Phase with beagle
module load java; \java -Xmx4g -jar /data/tunglab/tpv/Programs/beagle.08Jun17.d8b.jar gt=anubisSW_hwe.recode.vcf.gz out=anubisSW_beag ne=40000; tabix anubisSW_beag.vcf.gz

# Break to individual chromosomes
for f in `seq 1 20`; do vcftools --gzvcf anubisSW_beag.vcf.gz --chr $f --recode --out $f.anubisSW ; bgzip $f.anubisSW.recode.vcf; tabix $f.anubisSW.recode.vcf.gz; done
# Get individual-vcfs for each haplotype
for f in `seq 1 20`; do for g in `cat ../00_anu_SW.list`; do sed -e s/CHROM/$f/g run.01a.panubis1.sh | sed -e s/INDIV/$g/g > g.sh; sbatch g.sh; rm g.sh; done; done
# Run ldhelmet
for f in `seq 1 20`; do sed -e s/CHROM/$f/g run.02a.panubis1.sh > r.all.sh; sbatch --mem=65000 --nice r.all.sh; rm r.all.sh; done

