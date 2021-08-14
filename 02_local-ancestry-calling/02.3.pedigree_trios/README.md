## Accuracy of local ancestry calling using pedigree trios

## ARIELLE, IF WE UPLOAD MAJ RULE CALLS, MAKE SURE COLUMN ORDER MATCHES 4a SCRIPTS

This directory contains multiple scripts for evaluating the quality of local ancestry calls from LCLAE using known parent-offspring trios from the Amboseli baboons (Supplementary Methods 7.3).

#### Get pedigree inconsistences for the pedigree inconsistencies model (scripts labelled 1-3...).

In R, run `run.01.prep_for_pedigree_inconsistencies.R` which generates two R data files (`local_ancestry_pedigree_trios_maskedSNPRCref.Rd` and `local_ancestry_pedigree_trios_unmaskedWallref.Rd`) containing tracts, pedigree trio info, list of pedigree individuals, and genomic positions for evaluating the consistent of ancestry calls within pedigree trios for two sets of ancestry tracts (one generated using the SNPRC reference panel, one generated using the Wall et al. 2016 low coverage reference panel supplemented with high coverage Mikumi sequences). These R data files can then uploaded to a computing cluster for parallelization across chromosomes (which will make things run much faster).

In a directory on a computing cluster containing both R data files, run the R scripts `run.02.get_ancestry_calls_SNPRCref.R` and `run.02.get_ancestry_calls_Wallref.R` using the commands below in order to generate chromosome-specific scripts. These scripts (which only differ in the ancestry calls they are using) will generate chromosome-specific files containing the ancestry calls ancestry calls for each individual in the pedigree analysis at each of our focal positions (positions data frame created in `run.01.prep_for_pedigree_inconsistencies.R`):

```console 
# R must be loaded in your computing environment (e.g., module load R, activate a conda environment with R loaded, etc.)
for f in `cat autosomes.list`; do sed -e s/CHROMOSOME/$f/g run.02.get_ancestry_calls_SNPRCref.R > $f.sh; sbatch --mem=300 $f.sh; rm $f.sh; done
for f in `cat autosomes.list`; do sed -e s/CHROMOSOME/$f/g run.02.get_ancestry_calls_Wallref.R > $f.sh; sbatch --mem=300 $f.sh; rm $f.sh; done
# can check that all scripts ran by looking for "done" written in the output files
grep "done" slurm* | wc -l #40 = the number of total scripts we ran so everything ran to completion

# Concatenate the output of the chromosome-specific files into a single file for each set of ancestry calls (i.e., SNPRC or Wall et al. reference panels)
# SNPRC
# Remove the header of each file
for file in ancestry_calls_maskedSNPRCref*; do tail -n +2 "$file" >> "nohead.$file" ; done
# Grab the header from one of the original files with a header (all files should have the same header so it should not matter which chromosome file you choose)
head -1 ancestry_calls_maskedSNPRCref_pedigree_trios_35kbpos_chr1.txt >> header
# Concatenate the header and all of the header-less files into a single file
cat header nohead.* >> all.ancestry_calls_maskedSNPRCref_pedigree_trios_35kbpos.txt; rm *head*
# Wall
for file in ancestry_calls_unmaskedWallref*; do tail -n +2 "$file" >> "nohead.$file"; done
head -1 ancestry_calls_unmaskedWallref_pedigree_trios_35kbpos_chr1.txt >> header
cat header nohead.* >> all.ancestry_calls_unmaskedWallref_pedigree_trios_35kbpos.txt; rm *head*

# Delete chromsome-specific files as these results are now stored in either the all files generated above
rm ancestry*

# Check that we have the expected number of total lines per set of ancestry calls (should equal the total number of positions + header = 73975 + 1 = 73976)
wc -l all*
#73976 all.ancestry_calls_maskedSNPRCref_pedigree_trios_35kbpos.txt
#73976 all.ancestry_calls_unmaskedWallref_pedigree_trios_35kbpos.txt
```
In the same directory, run the R scripts `run.03.get_ped_inconsistencies_SNPRCref.R` and `run.03.get_ped_inconsistencies_Wallref.R` to identify inconsistent ancestry calls
```console 
sbatch --mem=1G  run.03.get_ped_inconsistencies_SNPRCref.R
sbatch --mem=1G  run.03.get_ped_inconsistencies_Wallref.R

# can check that all scripts ran by looking for "done" written in the output files
grep "done" slurm* | wc -l #2 = the number of total scripts we ran so everything ran to completion
```

#### Get covariates for the pedigree inconsistencies models (scripts labelled with a 4...)

For each genomic window, we would also like to get information on the number of (1) ancestry informative markers, (2) FST, and (3) recombination rate which we will include as covariates in our models of pedigree inconsistencies. 

**(1) ANCESTRY INFORMATIVE MARKERS**: run the R script `run.04a.AIM_count_PANEL.R` using the command below in order to generate chromosome-specific scripts. This script requires majority rule ancestry calls across Amboseli individuals, split into chromosome specific files. To generate these files, take the `amboseli.majrule.txt` output generated in Section 02.1 from the appropriate reference panel (e.g., masked SNPRC, unmasked Wall et al.), and split into chromosome-specific files by:
```console
awk '{OFS="\t"; print >> ("aims_maskedSNPRCref_" $1 ".txt")}' amb.majrule.maskedSNPRCref.txt # generates chromosome-specific files of ancestry calls across all Amboseli individuals 
awk '{OFS="\t"; print >> ("aims_unmaskedWallref_" $1 ".txt")}' amb.majrule.unmaskedWallref.txt # generates chromosome-specific files of ancestry calls across all Amboseli individuals 

for f in `seq 1 20`; do sed -e s/CHROMOSOME/$f/g run.04a.AIM_count_SNPRCref.R > $f.sh; sbatch --mem=30000 $f.sh; rm $f.sh; done
for f in `seq 1 20`; do sed -e s/CHROMOSOME/$f/g run.04a.AIM_count_Wallref.R > $f.sh; sbatch --mem=30000 $f.sh; rm $f.sh; done

# check that all scripts ran by looking for "done" written in the output files
grep "done" slurm* | wc -l #40 = the number of total scripts we ran so everything ran to completion

# concatenate the output of the chromosome-specific files into a single file for each set of ancestry calls (i.e., SNPRC or Wall et al. reference panels)
cat for_pedigree_trios_AIM_count_maskedSNPRCref_pedigree_trios_chr* >> all.for_pedigree_trios_AIM_count_maskedSNPRCref_pedigree_trios.txt
cat for_pedigree_trios_AIM_count_unmaskedWallref_pedigree_trios_chr* >> all.for_pedigree_trios_AIM_count_unmaskedWallref_pedigree_trios.txt

# Check that we have the expected number of total lines per set of ancestry calls (should equal the total number of positions = 73975)
wc -l all.for_pedigree_trios_AIM_count*
#73975 all.for_pedigree_trios_AIM_count_maskedSNPRCref_pedigree_trios.txt
#73975 all.for_pedigree_trios_AIM_count_unmaskedWallref_pedigree_trios.txt
  
# Remove chromosome-specific files
rm for_pedigree_trios_AIM_count*txt
```

**(2) FST**: use the scripts `run.04b.FST_SNPRCref.sh` and `run.04b.FST_Wallref.sh` to calculate FST using vcftools with 35 kb windows and 500 bp step size
```console
sbatch --mem=500 run.04b.FST_SNPRCref.sh
sbatch --mem=500 run.04b.FST_Wallref.sh
```

**(3) RECOMBINATION RATE**: run the R script from `run.04c.Recomb.R` which is almost identical to the code in `Section 05, r01`, but modified for the input for this analysis, the positions in this analysis (in the positions data frame), and the output (we want the lengths used).
```console
sbatch --mem=200 run.04c.Recomb.R
```

#### Results from the pedigree inconsistencies (script labelled with a 5...)

Finally, in a directory containing all files generated up to this point, run `run.05.pedigree_inconsistencies_results.R` which calculates, compares, and plots the proportion of ancestry state inconsistencies using ancestry calls from two reference panels. It also runs a model to evaluate what predicts ancestry state inconsistencies across sites. Finally, it evaluates the proportion of ancestry state inconsistencies per pedigree trio and whether the minimum coverage of individuals in the trio might contribute to it. 
