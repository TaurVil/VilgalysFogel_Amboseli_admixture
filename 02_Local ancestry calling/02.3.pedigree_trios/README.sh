# This directory contains multiple scripts for evaluating the quality of local ancestry calls from LCLAE using known parent-offspring trios from the Amboseli baboons (Supplementary Methods 7.3).

#############################################################################################################################
# Get pedigree inconsistences for the pedigree inconsistencies model (scripts labelled 1-3...).
#############################################################################################################################

# In R, run 1prep_for_pedigree_inconsistencies.R which generates two R data files (local_ancestry_pedigree_trios_maskedSNPRCref.Rd and local_ancestry_pedigree_trios_unmaskedWallref.Rd) containing tracts, pedigree trio info, list of pedigree individuals, and genomic positions for evaluating the consistent of ancestry calls within pedigree trios for two sets of ancestry tracts (one generated using the SNPRC reference panel, one generated using the Wall et al. 2016 Molecular Ecology low coverage reference panel). These R data files can then uploaded to a computing cluster for parallelization across chromosomes (which will make things run much faster).

# In a directory on a computing cluster containing both R data files, run the R scripts 2get_ancestry_calls_SNPRCref.R and 2get_ancestry_calls_Wallref.R using the commands below in order to generate chromosome-specific scripts. These scripts (which only differ in the ancestry calls they are using) will generate chromosome-specific files containing the ancestry calls ancestry calls for each individual in the pedigree analysis at each of our focal positions (positions data frame created in 1prep_for_pedigree_inconsistencies.R):
# R must be loaded in your computing environment (e.g., module load R, activate a conda environment with R loaded, etc.)
for f in `cat autosomes.list`; do sed -e s/CHROMOSOME/$f/g 2get_ancestry_calls_SNPRCref.R > $f.sh; sbatch --mem=300 $f.sh; rm $f.sh; done
for f in `cat autosomes.list`; do sed -e s/CHROMOSOME/$f/g 2get_ancestry_calls_Wallref.R > $f.sh; sbatch --mem=300 $f.sh; rm $f.sh; done
# can check that all scripts ran by looking for "done" written in the output files
grep "done" slurm* | wc -l #40 = the number of total scripts we ran so everything ran to completion

# Concatenate the output of the chromosome-specific files into a single file for each set of ancestry calls (i.e., SNPRC or Wall et al. reference panels)
# SNPRC
# Remove the header of each file
for file in ancestry_calls_maskedSNPRCref*
do
     tail -n +2 "$file" >> "nohead.$file" 
done

# Grab the header from one of the original files with a header (all files should have the same header so it should not matter which chromosome file you choose)
head -1 ancestry_calls_maskedSNPRCref_pedigree_trios_35kbpos_chr1.txt >> header
# Concatenate the header and all of the header-less files into a single file
cat header nohead.* >> all.ancestry_calls_maskedSNPRCref_pedigree_trios_35kbpos.txt
rm *head*

# Wall
for file in ancestry_calls_unmaskedWallref*
do
     tail -n +2 "$file" >> "nohead.$file" 
done

# Grab the header from one of the original files with a header (all files should have the same header so it should not matter which chromosome file you choose)
head -1 ancestry_calls_unmaskedWallref_pedigree_trios_35kbpos_chr1.txt >> header
# Concatenate the header and all of the header-less files into a single file
cat header nohead.* >> all.ancestry_calls_unmaskedWallref_pedigree_trios_35kbpos.txt
rm *head*

# Delete chromsome-specific files as these results are now stored in either the all files generated above
rm ancestry*

# Check that we have the expected number of total lines per set of ancestry calls (should equal the total number of positions + header = 73975 + 1 = 73976)
wc -l all*
#73976 all.ancestry_calls_maskedSNPRCref_pedigree_trios_35kbpos.txt
#73976 all.ancestry_calls_unmaskedWallref_pedigree_trios_35kbpos.txt

# In the same directory on a computing cluster, run the R scripts 3get_ped_inconsistencies_SNPRCref.R and 3get_ped_inconsistencies_Wallref.R to identify inconsistent ancestry calls
sbatch --mem=1G  3get_ped_inconsistencies_SNPRCref.R
sbatch --mem=1G  3get_ped_inconsistencies_Wallref.R

#############################################################################################################################
# Covariates for the pedigree inconsistencies models (scripts labelled with a 4...)
#############################################################################################################################

# For each genomic window, we would also like to get information on the number of ANCESTRY INFORMATIVE MARKERS, FST, and RECOMBINATION RATE which we will include as covariates in our models of pedigree inconsistencies
# For the number of ancestry informative markers, run the R script 4apos_AIM_count.sh using the command below in order to generate chromosome-specific scripts.
for f in `seq 1 20`; do sed -e s/CHROMOSOME/$f/g 4aAIM_count_SNPRCref.R > $f.sh; sbatch --mem=30000 $f.sh; rm $f.sh; done
for f in `seq 1 20`; do sed -e s/CHROMOSOME/$f/g 4aAIM_count_Wallref.R > $f.sh; sbatch --mem=30000 $f.sh; rm $f.sh; done

# Get first four columns (chrom, pos, start_window, end_window)
awk '{print $1,$2,$3,$4}' indiv.pos_AIM_count_1kbwin_HAP.txt >> tmp
for f in `ls indiv.pos_AIM_count_1kbwin_*`; do awk '{print $5}' $f >> tmp1.$f; echo $f; done
for f in `ls tmp1.indiv.pos_AIM_count_1kbwin_*`; do paste $f >> tmp1.$f; echo $f; done

# For FST, use the scripts 4bfst_SNPRCref.sh and 4bfst_Wallref.sh to calculate FST using vcftools with 35 kb windows and 500 bp step size
sbatch --mem=500 4bfst_SNPRCref.sh
sbatch --mem=500 4bfst_Wallref.sh

# For recombination rate, use the script XXXXTAURAS DIRECTORYXXX modified for the positions and genomic windows for this analysis (in the positions data frame).

#############################################################################################################################
# Results from the pedigree inconsistencies models (script labelled with a 5...)
#############################################################################################################################


# check Ns in chr7
# perl -ne 'chomp;if( />(.*)/){$head = $1; $i=0; next};@a=split("",$_); foreach(@a){$i++; if($_ eq "N" && $s ==0 ){$z=$i-1; print "$head\t$z"; $s =1}elsif($s==1 && $_ ne "N"){$j=$i-1;print "\t$j\n";$s=0}}' chr7.test.fa >> check_Ns
# SNP count
# for f in `seq 1 20`; do awk -v var=$f '{OFS="\t";if ($1==var) print }' SNP_count_35kbwin.$f.txt > SNP_count_35kbwin.$f.final.txt; done


# Finally, in R, run 3permuted_pvalue_calculation.R which calculates a permutation based p-value for each predictor variable based on the number of times that the absolute value of the effect size estimated from the permuted data sets was greater than the absolute value of the effect size estimated from the observed data set, across 1,000 permutations. 
