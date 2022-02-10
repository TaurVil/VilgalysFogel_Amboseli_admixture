## Genotype calling

#### Part 1: Scripts to process whole genome short-read Illumina sequence data (fastq -> BAM -> gVCF -> VCF).

For each sample, sequentially run WGSproc steps 1-7 to trim, map, and call genotypes per individual and then perform joint genotype calling. These jobs can be submitted together using run.baboonWGS_wynton_remotesubmit_20190814_single_highcov.sh. 

Joint genotype files are titled `baboon_genotypes.raw.CHROMOSOME.vcf.gz` and are available upon request or can be generated from fastq files uploaded to NCBI, but are not uploaded due to file size constraints for processed files. We are also happy to provide subsets of the dataset (i.e. just yellow or anubis baboons) upon request. Note that we performed joint genotyping with other samples that were subsequently removed before downstream analysis detailed in this directory.  

#### Part 2: Produce filtered files of anubis, yellow, and Amboseli baboon ancestry calls

First, we get genotype calls for putatively unadmixed yellow and anubis baboons. This requires unfiltered VCF files of baboon genotypes (`baboon_genotypes.raw.CHROMOSOME.vcf.gz`; called above and available upon request) and lists of unadmixed individuals (`00_anubis.list` and `00_yellow.list`, which can be found in Table S1). We hard filter these genotypes for sites called in >50% of both yellow and anubis baboons, as well as sites that pass hard filtering criteria (GATK best practices). Finally, we remove singleton and doubleton variants and export files of yellow and anubis genotype calls as `02.yel.$chrom.recode.vcf.gz` and `02.anu.$chrom.recode.vcf`. These files are combined into `yellow.vcf.gz` and `anubis.vcf.gz`, and merged into `refpanel.vcf.gz`. 

```console
## get calls for unadmixed individuals
sbatch --array=1-20 --mem=12G ./run.01.get_nonAmboseli_calls.sh

## merge vcfs for unadmixed yellow and anubis baboons
module load bcftools
bcftools concat ./chrom_vcfs/02.anu.*.vcf.gz -O z -o ./anubis.vcf.gz; tabix ./anubis.vcf.gz
bcftools concat ./chrom_vcfs/02.yel.*.vcf.gz -O z -o ./yellow.vcf.gz; tabix ./yellow.vcf.gz

## create merged file
bcftools merge ./anubis.vcf.gz ./yellow.vcf.gz -O z -o ./refpanel.vcf.gz

## get the mean depth per individual for each vcf
module load vcftools
vcftools --gzvcf ./anubis.vcf.gz --depth --out anubis
vcftools --gzvcf ./yellow.vcf.gz --depth --out yellow
```

Next we will get genotype calls for Amboseli individuals. This will require the unfiltered VCF file of baboon genotypes (will be made available on Zenodo as `baboon_genotypes.raw.CHROMOSOME.vcf.gz`) and a list of the Amboseli baboons (`00_amboseli.list`, which can be found in Table S1). We will retain biallelic variants which pass hard filtering criteria and are called in at least 20% of Amboseli individuals. We will not apply a minor allele frequency threshold within Amboseli, but rather focus on variants segregating between anubis and yellow baboons. These files are then combined into `amboseli.vcf.gz`. 

```console
sbatch --array=1-20 --mem=12G ./run.02.get_Amboseli_calls.sh

bcftools concat ./chrom_vcfs/02.*.amboseli.recode.vcf.gz -o ./amboseli.vcf.gz -O z; tabix ./amboseli.vcf.gz
```

Finally, we will create a merged call set of both Amboseli and unadmixed yellow and anubis individuals, including only samples that passed filtering criteria in both datasets. These files are saved by chromosome as `04.merged_shared.$index.vcf.gz` and merged into `merged_shared.vcf.gz`. 

```console
sbatch --array=1-20 --mem=12G ./run.03.merge.sh

bcftools concat ./chrom_vcfs/04.merged_shared.*.vcf.gz -O z -o ./merged_shared.vcf.gz
```

#### Part 3: Masking possibly introgressed anubis and yellow ancestry in non-Amboseli individuals

Our results revealed unexpected anubis ancestry in yellow baboons used to found the SNPRC baboon colony as well as low levels of apparently heterospecific ancestry in other unadmixed individuals, potentially consistent with background noise (see Section 4). We therefore sought to mask regions of the genome which contained putative introgressed ancestry in baboons outside of Amboseli, prior to estimate yellow and anubis baboon allele frequencies. To minimize the proportion of introgressed ancestry retained while also ensuring a high density of markers called in yellow and anubis individuals, we masked regions of the genome called as introgressed using LCLAE and those called IBD with heterospecific baboons, using IBDmix with at least 50% of possible individuals (see Section 4). 

Creating these files requires genotype calls for unadmixed individuals (`yellow.vcf.gz` and `anubis.vcf.gz`, generated above) and tracts of sequence to mask (`yes_intersect_50.bed`, generated in Section 04.2). These genotype calls will be available on Zenodo as `masked_yellow_and_anubis.vcf.gz`. 

```console
## we'll use tracts generated in section 04.2, specifically yes_intersect_50.bed. The title `yes_intersect_50` means (i) we are masking regions called introgressed in low coverage individuals using the results from LCLAE, (ii) we're using the intersection of LCLAE and IBDmix for high coverage individuals, and (iii) we're requiring IBD with 50% of source individuals (i.e. unadmixed yellow or anubis baboons) to include a region. 

## get tracts to mask per individual
for samp in `cat 00_refpanel.names`; do grep $samp yes_intersect_50.bed > ./tracts_to_mask/to_remove.$samp.yes_intersect_50.bed; done; echo $samp; done 
## replace empty files with a minimal bed file (just the first base, where no SNPs are called). This makes the files accessible to bedtools/vcftools down the road because they are in the correct format
for f in `find . -name '*.bed' -size 0`; do echo -e "chr1\t1\t2" > $f; done 

## get individual vcf files 
## calls a combined file of non-Amboseli yellow and anubis baboons created above: ./refpanel.vcf.gz
mkdir indiv_vcfs; for samp in `cat 00_refpanel.names`; do sed -e s/SAMPLE_NAME/$samp/g ./run.01.get_indiv_vcf.sh > ./r.$samp.sh; sbatch --mem=3G --nice ./r.$samp.sh; rm ./r.$samp.sh; done; done 

## merge individual masked vcfs back together
## filter for sites that still pass coverage and minor allele frequency filters
mkdir refpanel_vcfs; for vers in `echo yes_intersect_50`; do sed -e s/VERSION_NAME/$vers/g run.02.merge_filter_vcfs.sh > r.$vers.sh; sbatch --mem=8G r.$vers.sh; rm r.$vers.sh; done 
```

For calling local ancestry with LCLAE and Ancestry HMM, we merged genotype calls for Amboseli baboons with these masked genotype calls. This file (`amboseli_with_masked_refpanel.vcf.gz`) can be recreated from files on Zenodo using the code below. 

```console
## Amboseli genotype calls are available from Part 2 as ./amboseli.vcf.gz

## get list of sites in the masked, unadmixed baboon panel and that were called in amboseli baboons
## make sure both files have the same chromosome prefix (none or `chr`). If one doesn't, it can be removed from the other using `zcat MYVCF.vcf.gz  | sed -e 's/chr//g' | bgzip > ./OUT.vcf.gz`
vcftools --gzvcf ./masked/masked_yellow_and_anubis.vcf.gz --kept-sites --out ./masked/masked
vcftools --gzvcf ./amboseli.vcf.gz --kept-sites --out ./masked/amboseli

## make list of sites to keep (uniq -d deletes sites only seen in one file; grep -v 'CHROM' removes the header row)
sort ./masked/amboseli.kept.sites ./masked/masked.kept.sites | uniq -d | grep -v 'CHROM' > ./masked/kept.sites

## merge calls
## filter for sites called in both datasets
bcftools merge ./amboseli.vcf.gz ./masked/masked_yellow_and_anubis.vcf.gz -R masked_final/kept.sites -O z -o ./masked/amboseli_with_masked_refpanel.vcf.gz 
tabix ./masked/amboseli_with_masked_refpanel.vcf.gz

## split this into chromosome-specific vcf files
for chrom in {1..20}; do vcftools --gzvcf ./masked/amboseli_with_masked_refpanel.vcf.gz --chr $chrom  --recode --recode-INFO-all --out  amboseli_with_masked_refpanel.$chrom ; done
```
