## Local ancestry calling for Amboseli individuals using LCLAE

This folder contains instructions and scripts for generating local ancestry calls and tracts from genotype data using LCLAE (Wall et al. 2016 Molecular Ecology: https://doi.org/10.1111/mec.13684, https://github.com/jdwall02/LCLAE). You must download LCLAE for these scripts (available here: https://github.com/jdwall02/LCLAE)

FIrst, we create an input file of genotype calls from Amboseli baboons as well as yellow baboons and anubis baboons masked for putative introgression. This vcf (`amboseli_with_masked_refpanel.vcf.gz`) can be produced from `masked_yellow_and_anubis.vcf.gz`, which is available on Zenodo (doi://10.5281/zenodo.5199534; see section 01.3), and `amboseli.vcf.gz` (see section 01.2). We convert these genotype data into genotype likelihoods using LCLAE's `filtbaboon2b` function. In the script run.01.LCLAE_get_genolik.sh, XXX needs to be replaced with the number of individuals in the vcf file. 

```console 
# get list of individuals in the vcf file in the order that they appear in the vcf file
module load bcftools
bcftools query -l amboseli_with_masked_refpanel.vcf.gz >> 00_vcf_sample_order.masked.list 

# for each chromosome (1-20), format vcf for LCLAE and then use LCLAE's filtbaboon1b to convert genotype data into genotype likelihoods (run run.01.LCLAE_get_genolik.sh)
mkdir genotype_likelihoods_maskedref # make directory where we will store all of our genotype likelihood files
sbatch --array=1-20 --mem=100 run.01.LCLAE_get_genolik.sh
# note: in the script run.01.LCLAE_get_genolik.sh, on line 17 make sure to specify the total number of individuals in the vcf after calling the filtbaboon1b function (currently XXX is used as a placeholder for this number)

```

After generating genotype likelihood files, we define the reference individuals LCLAE will use (n=24 anubis and 7 yellow baboons used to found the SNPRC baboon colony), creating the files `ref_anubis.h` and `ref_yellow.h`. 

```console
# define our reference individuals
# for the two reference populations, separately get the order of each reference individual in the vcf file (e.g., are they individuals #1, #3, and #5?)
# here, we use the SNPRC yellow and anubis baboon founders as our reference individuals for yellow and anubis baboon ancestry, respectively. These lists are included in the ./DATA folder and are also listed in Table S1.
for i in `cat 00_anubis_SNPRC.list`; do grep -n $i 00_vcf_sample_order.masked.list >> ref_anubis; done
for i in `cat 00_yellow_SNPRC.list`; do grep -n $i 00_vcf_sample_order.masked.list >> ref_yellow; done
# confirm that we have the number of expected individuals in each reference panel (we will need these two numbers later)
wc -l ref_* # 24 ref_anubis; 7 ref_yellow
# grab the individuals' order numbers (we don't need their ids anymore)
sed -e s/\:.*//g ref_anubis >> ref_anubis2 
sed -e s/\:.*//g ref_yellow >> ref_yellow2
# sort the order numbers numerically
sort -n ref_anubis2 >> ref_anubis3
sort -n ref_yellow2 >> ref_yellow3
# output all order numbers so that they are all on a single line
echo $(cat ref_anubis3) >> ref_anubis.h
echo $(cat ref_yellow3) >> ref_yellow.h
# add the number of individuals per reference panel as the first entry in each line
sed -i -e 's/^/24 /' ref_anubis.h
sed -i -e 's/^/7 /' ref_yellow.h
```

Next, we call local ancestry using LCLAE's `filtbaboon2c` and `geno_lik2` functions, for each individual and each chromosome, using a sliding window approach. Again, in these scripts,the total number of individuals in the vcf needs to be specified after calling the filtbaboon2c function (currently XXX is used as a placeholder for this number). After running LCLAE for each individual-chromosome combination, we post-process the output files by: (1) adding the appropriate chromosome number as a column in each file, (2) merging all chromosomes for each individual into a single file per individual, (3) replacing the individual's sample number with their id in the file name, (4) adding the individual's id as a column in the file, (5) using majority rule for calling the ancestry state at each ancestry informative site, and (6) assembling ancestry tracts. The resulting ancestry tracts for all Amboseli individuals are available on Zenodo (doi://10.5281/zenodo.5199534) as `amboseli_LCLAE_tracts.txt`. 

```console 
# run LCLAE's filtbaboon2c and geno_lik2 to call local ancestry for each chromosome and individual
# NUMBER specifies the number of the individual in the VCF file (1, 2, 3, ..., XXX total number of individuals)
# the array/index specifies the chromosome for each iteration
# for example, to call local ancestry for all 20 chromosomes for the first 10 individuals in the vcf, we would run:
for f in `seq 1 10`; do sed -e s/NUMBER/$f/g run.02.LCLAE_get_ancestry.sh > r.$f.sh; sbatch --array=1-20 --mem=30000 r.$f.sh; rm r.$f.sh; done
# note: in the script run.02.LCLAE_get_ancestry.sh, on line 12, make sure to specify the total number of individuals in the vcf after calling the filtbaboon2c function (currently XXX is used as a placeholder for this number)
# returns an ancestry call at each ancestry informative SNP (minimum 20% allele frequency difference between the two reference populations), based on the surrounding genomic window (here, 35 kb)

# LCLAE sometimes fails, returning a "Segmentation Fault" error. If it works, the output (here, slurm-[job_id_#].out) will have a single output line with the sample number and the chromosome. 
# if not successful, there will be 4 rows 
# to identify the individual-chromosome files that failed, put those failed files into "redo2", and then re-run those individual-chromosome files.
# DO NOT RUN THIS UNTIL ALL JOBS HAVE FINISHED!
wc -l slurm* | grep ' 4 ' > redo; sed -i 's/[ ]*4 //g' redo; wc -l redo; for f in `cat redo`; do head -1 $f >> redo2; done; rm slurm*; wc -l redo2
tmp=`wc -l redo2 | awk '{print $1}'`
for i in `seq 1 $tmp`; do sed "${i}q;d" redo2 > tmp2; f=`awk '{print $1}' tmp2`; g=`awk '{print $2}' tmp2`; sed -e s/NUMBER/$f/g run.02b.rerun_LCLAE.sh | sed -e s/CHROMOSOME/$g/g > g.$f.$g.sh; sbatch g.$f.$g.sh; rm g.$f.$g.sh; done; rm redo*; rm tmp2
# each time, wait for all jobs to finish (!) and then repeat above step until all jobs successfully run (`wc -l redo` == 0).
# note: in the script run.02b.rerun_LCLAE.sh, on line 9, make sure to specify the total number of individuals in the vcf after calling the filtbaboon2c function (currently XXX is used as a placeholder for this number)
# note: another method to avoid re-running iterations which were already successful can be seen in Section 02.2. The segmentation fault error returns an output file with no rows, so running LCLAE within an if-then statement for the output file being 0 rows only runs iterations which failed. 

# if you get 0 redo (!), do the following as a second check to make sure all individual-chromosome files ran successfully
wc -l *35kb*txt > n_calls # get number of lines/calls for each individual-chromosome file
sed -i 's/^[ \t]*//' n_calls # fixes formatting of n_calls
grep '^0' n_calls > no_calls # grab files with 0 lines (that means something failed)
sed -i 's/\./\t/g' no_calls
wc -l no_calls # how many files have 0 lines - if 0, move on to post-processing
#0 no_calls
rm *n*calls
```

```console 
# add chromosome number as a column to each file, merge all chromosomes for each individual into a single file per individual, rename each individual's file with their id (from 00_vcf_sample_order.masked.list) as opposed to their number in the vcf file, and add the individual's id as a column in the file - note that there are several ways one could do these steps but below are the basic steps

# as above, "NUMBER" corresponds to the number of the individual in the vcf file (1, 2, 3, ..., XXX total number of individuals)

# for example, if we called local ancestry for the first 10 individuals in the vcf, we would then run:
for f in `seq 1 10`; do sed -e s/NUMBER/$f/g run.03.add_chrom_number.sh > g.$f.sh; sbatch --mem=16000 g.$f.sh; done # adds column with chromosome number to each file
rm g.*.sh
sbatch run.04.concatenate_files_per_indiv.sh # concatenates all chromosome files for each individual into a single file per individual
for h in `seq 1 10`; do tmp=`head -$h 00_vcf_sample_order.masked.list  | tail -1`; mv $h.35kb.d2.masked.SWref.txt $tmp.35kb.d2.masked.SWref.txt; done # change file name so it corresponds to the individual's actual id and not their number in the vcf file
sbatch run.05.add_indiv_id.sh # add column with the individual's id

# finally, for each ancestry informative site, use majority rule to assign the likely ancestry state using all sites within 35 kb of that site, and assemble ancestry tracts by combining contiguous sites that have the same ancestry assignment
# run R script run.06.majcalls_tracts.R for each individual
# if you generated ancestry calls for all individuals in the vcf and now want to use majority rule for ancestry assignment and tracts, do:
for f in `cat 00_vcf_sample_order.masked.list`; do cat run.06.majcalls_tracts.R  | sed -e s/INDIV/$f/g > g.$f.sh; sbatch --mem=45000 g.$f.sh; done
rm g.*.sh

# merge tracts for all amboseli individuals into a single file
touch amboseli_LCLAE_tracts.txt; for f in `ls *tracts.d2.*`; do sed '1d' $f >> amboseli_LCLAE_tracts.txt; done
# merge all majority rule calls for all individuals into a single file
touch amboseli.majrule.txt; for f in `ls *Maj*`; do sed '1d' $f >> amboseli.majrule.txt; done # these are only used in Section 02.3 for pedigree inconsistencies
```

## Visualizing local ancestry calls

We visualize local ancestry results using LCLAE in Figures 1A, 3B, and 3C. The code for each of these figures is provided here. 

* **Fig 1A**: visualization of local ancestry for the first 20 Mb of chromosome 1 for a random subsample of 100 Amboseli individuals. Calls data from `amboseli_LCLAE_tracts.txt`.
* **Fig 3B**: visualization of local ancestry genome-wide for two Amboseli individuals with different histories of admixture (e.g., recent vs. historically admixed). Calls data from `amboseli_LCLAE_tracts.txt`.
* **Fig 3C**: histogram of genome-wide anubis ancestry for the Amboseli baboons. Calls genome-wide anubis ancestry from Table S1 or the ids dataframe from `VilgalysFogel_main_data_file.250kb_windows.RData`. 
