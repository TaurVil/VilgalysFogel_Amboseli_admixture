# This directory contains multiple scripts for evaluating the quality of local ancestry calls from LCLAE using known parent-offspring trios from the Amboseli baboons (Supplementary Methods 7.3).

# In R, run 1prep_for_pedigree_inconsistencies.R which generates two R data files (local_ancestry_pedigree_trios_maskedSNPRCref.Rd and local_ancestry_pedigree_trios_unmaskedWallref.Rd) containing tracts, pedigree trio info, list of pedigree individuals, and genomic positions for evaluating the consistent of ancestry calls within pedigree trios for two sets of ancestry tracts (one generated using the SNPRC reference panel, one generated using the Wall et al. 2016 Molecular Ecology low coverage reference panel). These R data files can then uploaded to a computing cluster for parallelization across chromosomes (which will make things run much faster).

# In a directory on a computing cluster containing both R data files, run the R scripts 2get_ancestry_calls_SNPRCref.R and 2get_ancestry_calls_Wallref.R using the commands below in order to generate chromosome-specific scripts. These scripts (which only differ in the ancestry calls they are using) will generate chromosome-specific files containing the ancestry calls ancestry calls for each individual in the pedigree analysis at each of our focal positions (positions data frame created in 1prep_for_pedigree_inconsistencies.R):
# R must be loaded in your computing environment (e.g., module load R, activate a conda environment with R loaded, etc.)
for f in `cat autosomes.list`; do sed -e s/CHROMOSOME/$f/g 2get_ancestry_calls_SNPRCref.R > $f.sh; sbatch --mem=300 $f.sh; rm $f.sh; done
for f in `cat autosomes.list`; do sed -e s/CHROMOSOME/$f/g 2get_ancestry_calls_Wallref.R > $f.sh; sbatch --mem=300 $f.sh; rm $f.sh; done
# can check that all scripts ran by looking for "done" written in the output files
grep "done" slurm* | wc -l #40 = the number of total scripts we ran so everything ran to completion

# Concatenate the output of the chromosome-specific files into a single file for set of ancestry calls (i.e., SNPRC or Wall et al. reference panels)
for file in ancestry_calls_maskedSNPRCref*
do
     tail -n +2 "$file" >> "nohead.$file.txt" 
done

# Grab the header from one of the files
head -1 local_ancestry_calls_maskedfullref_pedigree_35kbpos_chr1.txt >> header
# Add the header to the concatenated file
cat header nohead.local_ancestry_calls_maskedfullref_pedigree_35kbpos_chr* >> all.local_ancestry_calls_maskedfullref_pedigree_35kbpos.txt
rm *head*

# we can also check that we have the expected number of total lines per set of ancestry calls (should equal the total number of positions = 73975)
wc -l all*


# get pedigree inconsistencies
sbatch --mem=1G run.02.ped_inconsistencies.sh

# get AIM count
for f in `seq 1 20`; do sed -e s/CHROMOSOME/$f/g 3pos_AIM_count.sh > $f.sh; sbatch --mem=30000 $f.sh; rm $f.sh; done



# get first four columns (chrom, pos, start_window, end_window)
awk '{print $1,$2,$3,$4}' indiv.pos_AIM_count_1kbwin_HAP.txt >> tmp
for f in `ls indiv.pos_AIM_count_1kbwin_*`; do awk '{print $5}' $f >> tmp1.$f; echo $f; done
for f in `ls tmp1.indiv.pos_AIM_count_1kbwin_*`; do paste $f >> tmp1.$f; echo $f; done

# check Ns in chr7
perl -ne 'chomp;if( />(.*)/){$head = $1; $i=0; next};@a=split("",$_); foreach(@a){$i++; if($_ eq "N" && $s ==0 ){$z=$i-1; print "$head\t$z"; $s =1}elsif($s==1 && $_ ne "N"){$j=$i-1;print "\t$j\n";$s=0}}' chr7.test.fa >> check_Ns

# SNP count

for f in `seq 1 20`; do awk -v var=$f '{OFS="\t";if ($1==var) print }' SNP_count_35kbwin.$f.txt > SNP_count_35kbwin.$f.final.txt; done




# Concatenate the output of the 50 jobs per array (20 permutations per job) into a single file for each model
cat groom_permuted_results_*.txt >> all.groom_permuted_results.txt
cat prox_permuted_results_*.txt >> all.prox_permuted_results.txt

# Delete individual files for each job in the array as these results are now stored in either the all.groom_permuted_results.txt or all.prox_permuted_results.txt files generated above
rm groom_permuted_results_*
rm prox_permuted_results_*

# Check that each dataframe with all of the permutation results has 1,000 lines corresponding to the 1,000 permutations run
wc -l all.*
#1000 all.groom_permuted_results.txt
#1000 all.prox_permuted_results.txt
#2000 total

# Finally, in R, run 3permuted_pvalue_calculation.R which calculates a permutation based p-value for each predictor variable based on the number of times that the absolute value of the effect size estimated from the permuted data sets was greater than the absolute value of the effect size estimated from the observed data set, across 1,000 permutations. 
