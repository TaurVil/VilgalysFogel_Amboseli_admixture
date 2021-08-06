### This folder contains instructions and scripts for generating local ancestry calls/tracts from genotype calls using LCLAE (Wall et al. 2016 Molecular Ecology: https://doi.org/10.1111/mec.13684, https://github.com/jdwall02/LCLAE ###

# Use vcfs (split by chromosome) masked for putative introgression available at XXXX
# Get list of individuals in the vcf file in the order that they appear in the vcf file (all chromosome vcf files have the same set of individuals so just grab the list of ordered individuals from chromosome 1)
module load bcftools
bcftools query -l analyzed_sites.masked_and_unmasked_refpanel.1.vcf.gz >> 00_vcf_sample_order.masked.list 

# For each chromosome (1-20), format vcf for LCLAE and then use LCLAE's filtbaboon1b to get genotype likelihoods (run 1LCLAE_get_genolik.sh)
mkdir genotype_likelihoods_maskedref # make directory where we will store all of our genotype likelihood files
sbatch --array=1-20 --mem=100 1LCLAE_get_genolik.sh
# Note: make sure the number after filtbaboon1b in run.04.get_genolik.sh is correct (the total number of individuals in the vcf)

# After generating genotype likelihood files, we'll call local ancestry using LCLAE's filtbaboon2c and geno_lik2 functions
# Before running filtbaboon3c and geno_lik2, we need to define our reference individuals
# For the two reference populations, separately get the order of each reference individual in the vcf file (e.g., are they individuals #1, #3, and #5?)
# Here, we'll use the SNPRC yellow and anubis baboon founders as our reference individuals for yellow and anubis baboon ancestry, respectively
for i in `cat 00_anubis_SNPRC.list`; do grep -n $i 00_vcf_sample_order.masked.list  >> ref_anubis; done
for i in `cat 00_yellow_SNPRC.list`; do grep -n $i 00_vcf_sample_order.masked.list >> ref_yellow; done

# Confirm that we have the number of expected individuals in each reference panel (we will need these two numbers later)
wc -l ref_*
#24 ref_anubis
#7 ref_yellow

# Grab the individual's order # (we do need their id anymore)
sed -e s/\:.*//g ref_anubis >> ref_anubis2 
sed -e s/\:.*//g ref_yellow >> ref_yellow2

# Sort the order #s numerically
sort -n ref_anubis2 >> ref_anubis3
sort -n ref_yellow2 >> ref_yellow3

# Output all numbers into a single line
echo $(cat ref_anubis3) >> ref_anubis.h
echo $(cat ref_yellow3) >> ref_yellow.h

# Add the number of individuals per reference panel as the first entry in each line
sed -i -e 's/^/24 /' ref_anubis.h
sed -i -e 's/^/7 /' ref_yellow.h

# Run LCLAE's filtbaboon3c and geno_lik2 to call local ancestry for each chromosome for each individual using a sliding window approach (here, 35 kb) (run 2LCLAE_get_ancestry.sh)
# "NUMBER" corresponds to the number of the individual in the vcf file (1, 2, 3, ..., XXX total number of individuals)
# The array/index specifies the chromosome
# For example, if we want to call local ancestry for all 20 chromosomes for the first 10 individuals in the vcf, we would run:
for f in `seq 1 10`; do sed -e s/NUMBER/$f/g 2LCLAE_get_ancestry.sh > r.$f.sh; sbatch --array=1-20 --mem=30000 r.$f.sh; done
rm r.$f.sh
# Returns an ancestry call at each ancestry informative SNP (minimum 20% allele frequency difference between the two reference populations), based on the surrounding genomic window (here, 35 kb). 

# LCLAE sometimes fails, returning a "Segmentation Fault" error. If it works, the output (here, slurm-[job_id_#].out) will have a single output line with the sample number and the chromosome. 
# If not successful, there will be 4 rows. 
# To identify the individual-chromosome files that failed, put those failed files into "redo2", and then re-run those individual-chromosome files.
# DO NOT RUN THIS UNTIL ALL JOBS HAVE FINISHED!
wc -l slurm* | grep ' 4 ' > redo; sed -i 's/[ ]*4 //g' redo; wc -l redo; for f in `cat redo`; do head -1 $f >> redo2; done; rm slurm*; wc -l redo2
tmp=`wc -l redo2 | awk '{print $1}'`
for i in `seq 1 $tmp`; do sed "${i}q;d" redo2 > tmp2; f=`awk '{print $1}' tmp2`; g=`awk '{print $2}' tmp2`; sed -e s/NUMBER/$f/g rerun.2LCLAE_get_ancestry.sh | sed -e s/CHROMOSOME/$g/g > g.$f.$g.sh; sbatch g.$f.$g.sh; rm g.$f.$g.sh; done; 
rm redo*; rm tmp2
# Each time, wait for all jobs to finish (!) and then repeat above step until all jobs successfully run (`wc -l redo` == 0).

# If you get 0 redo, do this as a second check to make sure all individual-chromosome files ran successfully
wc -l *35kb*txt > n_calls # get number of calls per chromosome
sed -i 's/^[ \t]*//' n_calls # fixes formatting of n_calls
grep '^0' n_calls > no_calls
sed -i 's/\./\t/g' no_calls
wc -l no_calls 
#0 no_calls
rm *n*calls

# Add chromosome number to each file, merge all chromosomes for each individual, and replace the individual's sample number with the individual's original file name (from 00_vcf_sample_order.masked.list) - note that there are several ways one could do these steps but below are the basic steps
# For example, if we want to call local ancestry for all 20 chromosomes for the first 10 individuals in the vcf, we would run:
# As above, "NUMBER" corresponds to the number of the individual in the vcf file (1, 2, 3, ..., XXX total number of individuals)
# For example, if we called local ancestry for the first 10 individuals in the vcf, we would run:
for f in `seq 1 10`; do sed -e s/NUMBER/$f/g 3add_chrom_number.sh > g.$f.sh; sbatch --mem=16000 g.$f.sh; rm g.$f.sh; done # adds column with chromosome number to each file
sbatch 4concatenate_files_per_indiv.sh # concatenates all chromosome files for each individual into a single file per individual
for h in `seq 1 10`; do tmp=`head -$h 00_vcf_sample_order.masked.list  | tail -1`; mv $h.35kb.d2.masked.SWref.txt $tmp.35kb.d2.masked.SWref.txt; done # change file name so it corresponds to the actual individual's id and not their number in the vcf file
sbatch 5add_indiv_id.sh # add column with the individual ID

# Finally, for each ancestry informative SNP, use majority rule to assign the likely ancestry state using all SNPs within 35kb of that site, and assemble ancestry tracts by combining contiguous SNPs that have the same ancestry assignment
# Run R script 7majcalls_tracts.R for each individual
# If you generated ancestry calls for all individuals in the vcf and now want to use majority rule for ancestry assignment and tracts, do:
for f in `cat 00_vcf_sample_order.masked.list`; do cat 7majcalls_tracts.R  | sed -e s/INDIV/$f/g > g.$f.sh; sbatch --mem=45000 g.$f.sh; done



## Merge all calls into a single file
touch all.1kb.tracts.fullref.txt; for f in `ls *1kb.d2.*`; do sed '1d' $f >> all.1kb.tracts.fullref.txt; done

# for fig 1,
touch all.majrule.maskedSWref_and_unmaskedref.txt; for f in `ls *Maj*`; do sed '1d' $f >>  all.majrule.maskedSWref_and_unmaskedref.txt; done

for f in `cat 00_amboseli.list`; do mv $f*.d2*txt ancestry_usingfullref/amboseli_indiv/; done
for f in `cat 00_allref.list`; do mv $f*d2*txt ancestry_usingfullref/refpanel_indiv/; done

# In amboseli_indiv directory, merge majority rule call and tract files
touch amb.tracts.txt; for f in `ls *tracts.d2.*`; do sed '1d' $f >> amb.tracts.txt; done
touch amb.majrule.txt; for f in `ls *Maj*`; do sed '1d' $f >> amb.majrule.txt; done

## Clean up step: get rid of genotype likelihoods, r.06_*, calls by chromosome, and any other intermediate files. ## clean up the tmp2.INDIV.txt files, leaving just the majority rule files and the tracts. 
