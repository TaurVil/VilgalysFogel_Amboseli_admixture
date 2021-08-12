## Genotype calling

#### Part 1: Scripts to process whole genome short-read Illumina sequence data (fastq -> BAM -> gVCF).

For each sample, sequentially run WGSproc steps 1-7 to trim, map, and call genotypes per individual. These jobs can be submitted together using baboonWGS_wynton_remotesubmit_20190814_single_highcov.sh. 

The resulting genotyping calls are available on Zenodo as XXX. 

Note that we performed joint genotyping with other samples that were subsequently removed before downstream analysis detailed in this directory. We don’t anticipate that the inclusion of other samples during the joint genotyping step affects our results given our subsequent filtering and processing before analysis. For any questions regarding this content, please contact Jacqueline Robinson (https://jarobinsonresearch.com). 

#### Part 2: Filter merged genotype calls to produce files of anubis and yellow baboon genotype calls, as well as a callset merged with the Amboseli data. 

First we'll get genotype calls for putatively unadmixed yellow and anubis baboons. This will require the unfiltered VCF file of baboon genotypes (made available on Zenodo as baboon_genotypes.raw.CHROMOSOME.vcf.gz) and files of unadmixed individuals (`00_anubis.list` and `00_yellow.list`, which can be found in Table S1). We'll hard filter these genotypes for sites called in >50% of both yellow and anubis baboons, as well as sites that pass hard filtering criteria. Finally, we'll remove singleton variants and export files of yellow and anubis genotype calls as `02.yel.$chrom.recode.vcf.gz` and `02.anu.$chrom.recode.vcf`. These files will then be combined into `yellow.vcf.gz` and `anubis.vcf.gz`, which are available on Zenodo. 

```console
## Get calls for unadmixed individuals
sbatch --array=1-20 --mem=12G ./r01.get_nonAmboseli_calls.sh

## merge vcfs for yellow and anubis baboons
module load bcftools
bcftools concat ./chrom_vcfs/02.anu.*.vcf.gz -O z -o ./anubis.vcf.gz; tabix ./anubis.vcf.gz
bcftools concat ./chrom_vcfs/02.yel.*.vcf.gz -O z -o ./yellow.vcf.gz; tabix ./yellow.vcf.gz

## create merged file
bcftools merge ./anubis.vcf.gz ./yellow.vcf.gz -O z -o ./refpanel.vcf.gz

## get the depth of coverage per SNP for each file
module load vcftools
vcftools --gzvcf ./anubis.vcf.gz --depth --out anubis
vcftools --gzvcf ./yellow.vcf.gz --depth --out yellow
```

Next we'll get genotype calls for Amboseli individuals. This will require the unfiltered VCF file of baboon genotypes (made available on Zenodo as baboon_genotypes.raw.CHROMOSOME.vcf.gz) and a file of Amboseli baboons to include (`00_amboseli.list`, which can be found in Table S1). We'll retain biallelic variants which pass hard filtering criteria and are called in at least 20% of Amboseli individuals. We will not apply a minor allele frequency threshold within Amboseli, but rather focus on variants segregating between anubis and yellow baboons. These files are then combined into `amboseli.vcf.gz`, which is available on Zenodo. 

```console
sbatch --array=1-20 --mem=12G ./r02.get_Amboseli_calls.sh

bcftools concat ./chrom_vcfs/02.*.amboseli.recode.vcf.gz -o ./amboseli.vcf.gz -O z; tabix ./amboseli.vcf.gz
```

Finally, we'll create a merged call set of both Amboseli and unadmixed individuals, including only samples that passed filtering criteria in both datasets. These files are saved by chromosome as `04.merged_shared.$index.vcf.gz` and merged into `merged_shared.vcf.gz`. 

```console
sbatch --array=1-20 --mem=12G ./r03.merge.sh

bcftools concat ./chrom_vcfs/04.merged_shared.*.vcf.gz -O z -o ./merged_shared.vcf.gz
```

#### Part 3: Masking anubis and yellow ancestry in non-Amboseli individuals

Our results revealed unexpected anubis ancestry in yellow baboons used to found the SNPRC as well as low levels of heterospecific ancestry in other unadmixed individuals, perhaps consistent with background noise (see Section 4). We therefore sought to mask regions of the genome which contained putative heterospecific ancestry from the individuals we would use as unadmixed reference panels when calling local ancestry and to estimate yellow and anubis baboon allele frequencies. To minimize the proportion of anubis ancestry maintained while also retaining a high density of markers called in sufficient yellow and anubis individuals, we masked tracts of the genome called as heterospecific using LCLAE and those IBD using at least 50% of possible source individuals (see Section 4). 

Creating these files will require genotype calls for unadmixed individuals (`yellow.vcf.gz` and `anubis.vcf.gz`, available on Zenodo) and tracts of sequence to mask (`yes_intersect_50.bed`, available in Section 04.2). These genotype calls are available on Zenodo titled `masked_yellow_and_anubis.vcf.gz`. 

```console
## we'll use tracts generated in section 04.2, specifically yes_intersect_50.bed. The title `yes_intersect_50` means (i) we are masking regions called heterospecific in low coverage individuals using the results from LCLAE, (ii) we're using the intersect of LCLAE and IBDmix for high coverage individuals, and (iii) we're requiring IBD with 50% of source individuals to include a region. 

## get tracts to mask per individual
for samp in `cat 00_refpanel.names`; do grep $samp yes_intersect_50.bed > ./tracts_to_mask/to_remove.$samp.yes_intersect_50.bed; done; echo $samp; done 
## Replace empty files with a minimal variant (just excluding the first BP, where no SNPs are called). This makes the files accessible to bedtools/vcftools down the road because they are in the correct format
for f in `find . -name '*.bed' -size 0`; do echo -e "chr1\t1\t2" > $f; done 

## get individual vcf files 
## Calls a combined reference panel of unadmixed individuals created above: ./refpanel.vcf.gz
mkdir indiv_vcfs; for samp in `cat 00_refpanel.names`; do sed -e s/SAMPLE_NAME/$samp/g ./run.01.get_indiv_vcf.sh > ./r.$samp.sh; sbatch --mem=3G --nice ./r.$samp.sh; rm ./r.$samp.sh; done; done 

## merge individual masked vcfs back together
## filter for sites that still pass coverage and minor allele frequency filters
mkdir refpanel_vcfs; for vers in `echo yes_intersect_50`; do sed -e s/VERSION_NAME/$vers/g run.02.merge_filter_vcfs.sh > r.$vers.sh; sbatch --mem=8G r.$vers.sh; rm r.$vers.sh; done 
```

For calling local ancestry with LCLAE and Ancestry HMM, we merged genotype calls for Amboseli baboons with these masked genotype calls. This file (`amboseli_with_masked_refpanel.vcf.gz`) can be recreated from files on Zenodo using the code below. 

```console
## Amboseli genotype calls are available from Part 2 as ./amboseli.vcf.gz

## Get sites in the refpanel and amboseli 
## make sure both files have the same chromosome prefix (none or `chr`). If one doesn't, it can be removed from the other using `zcat MYVCF.vcf.gz  | sed -e 's/chr//g' | bgzip > ./OUT.vcf.gz`
vcftools --gzvcf ./masked/masked_yellow_and_anubis.vcf.gz --kept-sites --out ./masked/masked
vcftools --gzvcf ./amboseli.vcf.gz --kept-sites --out ./masked/amboseli

## make list of sites to keep (uniq -d deletes sites only seen in one file; grep -v 'CHROM' removes the header row)
sort ./masked/amboseli.kept.sites ./masked/masked.kept.sites | uniq -d | grep -v 'CHROM' > ./masked/kept.sites

## merge calls
## filter for sites called in each dataset 
bcftools merge ./amboseli.vcf.gz ./masked/masked_yellow_and_anubis.vcf.gz -R masked_final/kept.sites -O z -o ./masked/amboseli_with_masked_refpanel.vcf.gz 
tabix ./masked/amboseli_with_masked_refpanel.vcf.gz

## split this into a vcf file per chromosome
for chrom in {1..20}; do vcftools --gzvcf ./masked/amboseli_with_masked_refpanel.vcf.gz --chr $chrom  --recode --recode-INFO-all --out  amboseli_with_masked_refpanel.$chrom ; done
```
