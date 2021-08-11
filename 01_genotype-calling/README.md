## Genotype calling

#### Part 1: Scripts to process whole genome short-read Illumina sequence data (fastq -> BAM -> gVCF).

For each sample, sequentially run WGSproc steps 1-7 to trim, map, and call genotypes per individual. These jobs can be submitted together using baboonWGS_wynton_remotesubmit_20190814_single_highcov.sh. 

The resulting genotyping calls are available on Zenodo as XXX. 

Note that we performed joint genotyping with other samples that were subsequently removed before downstream analysis detailed in this directory. We don’t anticipate that the inclusion of other samples during the joint genotyping step affects our results given our subsequent filtering and processing before analysis. For any questions regarding this content, please contact Jacqueline Robinson (https://jarobinsonresearch.com). 

#### Part 2: Filter merged genotype calls to produce files of anubis and yellow baboon genotype calls, as well as a callset merged with the Amboseli data. 

First we'll get genotype calls for putatively unadmixed yellow and anubis baboons. This will require the unfiltered VCF file of baboon genotypes (made available on Zenodo as XXX) and files of unadmixed individuals (`00_anubis.list` and `00_yellow.list`, which can be found in Table S1). We'll hard filter these genotypes for sites called in >50% of both yellow and anubis baboons, as well as sites that pass hard filtering criteria. Finally, we'll remove singleton variants and export files of yellow and anubis genotype calls as `bgzip 02.yel.$chrom.recode.vcf.gz` and `02.anu.$chrom.recode.vcf`. These files will then be combined into `yellow.vcf.gz` and `anubis.vcf.gz`. 

```console
## Get calls for unadmixed individuals
sbatch --array=1-20 --mem=12G ./r01.get_nonAmboseli_calls.sh

## merge vcfs for yellow and anubis baboons
module load bcftools
bcftools concat ./chrom_vcfs/02.anu.*.vcf.gz -O z -o ./anubis.vcf.gz
bcftools concat ./chrom_vcfs/02.yel.*.vcf.gz -O z -o ./yellow.vcf.gz

## get the depth of coverage per SNP for each file
module load vcftools
vcftools --gzvcf ./anubis.vcf.gz --depth --out anubis
vcftools --gzvcf ./yellow.vcf.gz --depth --out yellow

```

Next we'll get genotype calls for Amboseli individuals. This will require the unfiltered VCF file of baboon genotypes (made available on Zenodo as XXX) and a file of Amboseli baboons to include (`00_amboseli.list`, which can be found in Table S1). We'll retain biallelic variants which pass hard filtering criteria and are called in at least 20% of Amboseli individuals. We will not apply a minor allele frequency threshold within Amboseli, but rather focus on variants segregating between anubis and yellow baboons. 

```console
sbatch --array=1-20 --mem=12G ./r02.get_Amboseli_calls.sh

```

Finally, we'll create a merged call set of both Amboseli and unadmixed individuals, including only samples that passed filtering criteria in both datasets. These files are saved as `04.merged_shared.$index.vcf.gz`. 

```console
sbatch --array=1-20 --mem=12G ./r03.merge.sh
```