## Setup

# get yellow and anubis allele frequencies using SNPRC baboon colony founder yellow and anubis individuals (the same reference panel used for LCLAE)
# the vcf file used for these calls come from Section 01, and the process for getting allele frequencies is described in Section 03. We'll use the "frq" tables for chromosomes 17-20 here.
# these are saved in a_priori_genotypes, for use in the script run.02.create_sample_vcf.R

## get reduced genome for baboon chromsomes 17-20. This is the portion of the genome we will simulate ancestry and sequencing data for. 
mkdir partial_genome; module load samtools
for f in `cat 00chroms`; do samtools faidx ~/genomes/panubis1/Panubis1.0.fa $f > ./partial_genome/$f.fa; echo $f; done
cat ./partial_genome/chr*.fa > ./partial_genome/Reduced_Genome.fa
bowtie2-build ./partial_genome/Reduced_Genome.fa reduced_genome

## get header to be added to the vcf files after genotypes are simulated. This is added in run.03.get_fasta.sh, and can come from any vcf file. 
zcat ~/panubis1_genotypes/calls_merged/04.merged_shared.1.vcf.gz | grep '^#' | sed '$ d' > full_header.vcf 

## create last line of vcf header to which genotype information will be appended in run.02.create_sample_vcf.R
zcat ./calls_unadmixed/02.yel.20.recode.vcf.gz | grep '^#' | tail -1 > vcf_example_header
# manually edit to create 1 name only ("SAMPLE")

