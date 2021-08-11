
## DATES (Narasimhan et al. 2019, https://github.com/priyamoorjani/DATES) to estimate the timing of admixture in high coverage Amboseli baboons. 

Briefly, DATES (Distribution of Ancestry Tracts of Evolutionary Signals) is a method to estimate the time of admixture based on the decay in genetic ancestry between locations. It takes as inputs a genotype file (eigenstrat format), population to which each individual belongs, and a configuration file. It returns estimates for how many generations prior to the present admixture occurred. Dates relies upon high coverage genomes, so we focus on the 9 Amboseli individuals for whom we have >10x sequencing coverage.  


Because admixed individuals and missing genotype data can affect DATES, we used the high coverage Mikumi yellow baboons and SW anubis baboon founders to call DATES. Unlike with LCLAE, we were unable to use SW yellow founders because the missing genotypes from masking anubis ancestry would negatively affect DATES. We will get these genotypes from the combined vcf file, then filter for biallelic SNPs and convert to the ancestry map format.  

```console 
## Grab the right samples 
module load vcftools; module load plink
vcftools --keep 00_anu.list --keep 00_yel.list --vcf ../ALL_GENOTYPES.vcf.gz --recode --out tmp.ref.for_plink
plink --vcf tmp.ref.for_plink.recode.vcf --snps-only --maf 0.05 --recode --out plink.ref 
## Conversion from plink to ANCESTRYMAP format uses convertf: https://github.com/argriffing/eigensoft/tree/master/CONVERTF
~/Programs/EIG-6.1.4/bin/convertf -p par.PED.ANCESTRYMAP
## following this, the population column needs to be manually edited for yellow vs anubis baboons, which is saved as "ref.v2.ind"

## the snp file needs to be reformated to add chromosome positions in cM as well as the physical positions. 
## lengths of each chromosome in cM come from Cox et al. 2006 (object `d` in the R script)
module load R; ./r01.get_genetic_positions.R

```

We'll also pull the high coverage Amboseli genotypes and format those as well. 

```console
vcftools --keep 00_amboseli.list --vcf ../ALL_GENOTYPES.vcf.gz --recode --out tmp.ambo.for_plink
plink --vcf tmp.ambo.for_plink.recode.vcf --snps-only --recode --out plink.ambo 
~/Programs/EIG-6.1.4/bin/convertf -p par.PED.ANCESTRYMAP2

sed -e "s/ref/ambo/g" ./r01.get_genetic_positions.R > ./r01b.get_genetic_positions_ambo.R
./r01b.get_genetic_positions_ambo.R; rm r01b.get_genetic_positions_ambo.R

```

DATES provides a function to merge ANCESTRYMAP files, after which we need to do some reformatting to make sure SNPs are included in all sampels and sequential. 

```console
~/Programs/DATES-master/example/mergeit -p par.mergeit > mergeit.log 
./r02.fix_merged_file.R
~/Programs/EIG-6.1.4/bin/convertf -p par.convertf.order
```

Finally, we will run DATES itself to estimate the timing of admixture. 

```console
## Load relevant packages
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/Programs/gsl-2.3/lib:~/Programs/fftw-3.3.3/lib
export LD_LIBRARY_PATH; export PATH=$PATH:$LD_LIBRARY_PATH
export PATH=$PATH:~/Programs/DATES-master/src/bin/
module load OpenBLAS/0.2.20-gcb01; module load glibc/2.14-gcb01; module load gnuplot 

~/Programs/DATES-master/src/bin/dates -p par.dates > log.dates 
```

DATES outputs one folder per individual. Results are summarized in `DATES_results.txt` and plotted for Fig. S6 using `figureS6.R`

```console
./r03.figureS6.R
```