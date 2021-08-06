# Get list of individuals in the vcf file in the order that they appear in the vcf file (all chromosome vcf files have 509 individuals so just grab the list of ordered individuals from chromosome 1)
module load bcftools
bcftools query -l 04.merged_shared.1.vcf.gz >> 04_vcf_sample_order.list
bcftools query -l analyzed_sites.masked_and_unmasked_refpanel.1.vcf.gz >> 04_vcf_sample_order.masked_and_unmasked_ref.list

## Format vcf files for LCLAE from Wall et al. 2016 (https://github.com/jdwall02/LCLAE)
## The order of file names in the vcf file are in 04_vcf_sample_order.list

## Split vcf into chromosomes and get genotype likelihoods.
## There are 20 autosomes, an X, and a Y so should get 22 genotype likelihood files. Note that there's no mitochondrial genome in this assembly.
## Right now, only do autosomes
sbatch --array=1-20 --mem=100 run.03.LCLAE.get_genolik.sh
# Note: make sure the number after filtbaboon1b in run.04.get_genolik.sh is correct (the total number of individuals in the vcf). If not, everything gets screwed up. 
# Right now, this is 508 which is correct given our current dataset 

# for Figure 1, 
sbatch --array=1-20 --mem=100 run.03.LCLAE.get_genolik_3.sh 

# using the entire ref panel
# remove pcyn_16098 from 00_yel.list since they are a duplicate of Mik_07 according to Tauras
grep -v "pcyn_16098" 00_yel.list.nodup >> 00_yel.list.nodup.no16098
# get the order of each reference individual in the vcf file
 for i in `cat 00_anu.list`; do grep -n $i 04_vcf_sample_order.list >> fullref_anubis; done
 for i in `cat 00_yel.list.nodup.no16098`; do grep -n $i 04_vcf_sample_order.list >> fullref_yellow; done

# how many individuals in each ref panel? need this for later
wc -l fullref_*
#41 fullref_anubis
#22 fullref_yellow

# just grab the order # and not the individual name
sed -e s/\:.*//g fullref_anubis >> fullref_anubis2 
sed -e s/\:.*//g fullref_yellow >> fullref_yellow2

# sort numerically
sort -n fullref_anubis2 >> fullref_anubis3
sort -n fullref_yellow2 >> fullref_yellow3

# put output all into one line
echo $(cat fullref_anubis3) >> fullref_anubis.h
echo $(cat fullref_yellow3) >> fullref_yellow.h

# add number of individuals per reference panel as first entry 
 sed -i -e 's/^/41 /' fullref_anubis.h
 sed -i -e 's/^/22 /' fullref_yellow.h 

## using the wall ref panel
for i in `cat 00_anu_wall.list`; do grep -n $i 04_vcf_sample_order.list >> wallref_anubis; done
for i in `cat 00_yel_wall.list.nodup`; do grep -n $i 04_vcf_sample_order.list >> wallref_yellow; done 

wc -l wallref_*
#13 wallref_anubis
#15 wallref_yellow

sed -e s/\:.*//g wallref_anubis >> wallref_anubis2
sed -e s/\:.*//g wallref_yellow >> wallref_yellow2

# sort numerically
sort -n wallref_anubis2 >> wallref_anubis3
sort -n wallref_yellow2 >>  wallref_yellow3

echo $(cat wallref_anubis3) >> wallref_anubis.h
echo $(cat wallref_yellow3) >> wallref_yellow.h

sed -i -e 's/^/13 /' wallref_anubis.h
sed -i -e 's/^/15 /' wallref_yellow.h 

## using the SW ref panel
for i in `cat 00_anu_SW.list`; do grep -n $i 04_vcf_sample_order.list >> SWref_anubis; done
for i in `cat 00_yel_SW.list`; do grep -n $i 04_vcf_sample_order.list >> SWref_yellow; done 

wc -l SWref_*
#24 SWref_anubis
# 7 SWref_yellow

sed -e s/\:.*//g SWref_anubis >> SWref_anubis2
sed -e s/\:.*//g SWref_yellow >> SWref_yellow2

# sort numerically
sort -n SWref_anubis2 >> SWref_anubis3
sort -n SWref_yellow2 >>  SWref_yellow3

echo $(cat SWref_anubis3) >> SWref_anubis.h
echo $(cat SWref_yellow3) >> SWref_yellow.h

sed -i -e 's/^/24 /' SWref_anubis.h
sed -i -e 's/^/7 /' SWref_yellow.h 

## using the SW ref panel for figure 1
# masked SW individuals will be our reference but we're interested in plotting the ancestry of the calls of the unmasked ref individuals
sed -e s/^/masked_\/g 00_anu_SW.list >> 00_anu_SW.list2
sed -e s/^/masked_\/g 00_yel_SW.list >> 00_yel_SW.list2

for i in `cat 00_anu_SW.list2`; do grep -w -n $i 04_vcf_sample_order.masked_and_unmasked_ref.list >> SWref_anubis; done
for i in `cat 00_yel_SW.list2`; do grep -w -n $i 04_vcf_sample_order.masked_and_unmasked_ref.list >> SWref_yellow; done

wc -l SWref_*
#24 SWref_anubis
# 7 SWref_yellow

sed -e s/\:.*//g SWref_anubis >> SWref_anubis2
sed -e s/\:.*//g SWref_yellow >> SWref_yellow2

# sort numerically
sort -n SWref_anubis2 >> SWref_anubis3
sort -n SWref_yellow2 >>  SWref_yellow3

echo $(cat SWref_anubis3) >> maskedSWref_anubis_2.h
echo $(cat SWref_yellow3) >> maskedSWref_yellow_2.h

sed -i -e 's/^/24 /' maskedSWref_anubis_2.h
sed -i -e 's/^/7 /' maskedSWref_yellow_2.hh 

## Use LCLAE to call local ancestry for each individual using a sliding window (35 kb) approach
## Individual should correspond to the number of the individual in the vcf file (1, 2, 3, ..., 508)
for f in `seq 1 508`; do sed -e s/NUMBER/$f/g run.04a.get_ancestry_calls_masked_fullref.sh > r.04_$f.sh; sbatch --array=1-20%1 --mem=30000 r.04_$f.sh; done
rm r.04_*.sh
# Returns an ancestry call at each ancestry informative SNP, based on the surrounding 35kb window. 

# fig 1, use:
for f in `seq 1 126`; do sed -e s/NUMBER/$f/g run.04a.get_ancestry_calls_masked_SWref_and_unmasked_ref.sh > r.04_$f.sh; sbatch --array=1-20%3 --mem=30000 r.04_$f.sh; done
rm r.04_*.sh

## LCLAE sometimes fails, returning a "Segmentation Fault" error. If it works, the output will have a single output line with the sample number and the chromosome. 
## If not successful, there will be 4 rows. 
## We must identify the individual-chromosome files that failed, put those failed files into "redo2", and then re-run those individual-chromosome files.

## DO NOT RUN THIS UNTIL ALL JOBS HAVE FINISHED!!!
wc -l slurm* | grep ' 4 ' > redo; sed -i 's/[ ]*4 //g' redo; wc -l redo; for f in `cat redo`; do head -1 $f >> redo2; done; rm slurm*; wc -l redo2
tmp=`wc -l redo2 | awk '{print $1}'`
for i in `seq 1 $tmp`; do sed "${i}q;d" redo2 > tmp2; f=`awk '{print $1}' tmp2`; g=`awk '{print $2}' tmp2`; sed -e s/NUMBER/$f/g run.04b.rerun.sh | sed -e s/CHROMOSOME/$g/g > g.$f.$g.sh; sbatch g.$f.$g.sh; rm g.$f.$g.sh; done; 
rm redo*; rm tmp2
## Each time, wait for all jobs to finish (!) and then repeat above step until all jobs successfully run (`wc -l redo` == 0).

# arielle, can we remove this later - it's for you!
for k in `seq 1 508`; do for file in g.$k.*.sh; do sbatch $file; rm $file; done; done

# for fig 1,
wc -l slurm* | grep ' 4 ' > redo; sed -i 's/[ ]*4 //g' redo; wc -l redo; for f in `cat redo`; do head -1 $f >> redo2; done; rm slurm*; wc -l redo2
tmp=`wc -l redo2 | awk '{print $1}'`
for i in `seq 1 $tmp`; do sed "${i}q;d" redo2 > tmp2; f=`awk '{print $1}' tmp2`; g=`awk '{print $2}' tmp2`; sed -e s/NUMBER/$f/g run.04b.rerun_2.sh | sed -e s/CHROMOSOME/$g/g > g.$f.$g.sh; sbatch g.$f.$g.sh; rm g.$f.$g.sh; done; 
rm redo*; rm tmp2
## Each time, wait for all jobs to finish (!) and then repeat above step until all jobs successfully run (`wc -l redo` == 0).
for k in `seq 1 126`; do for file in g.$k.*.sh; do sbatch $file; rm $file;

# If you get 0 redo, do this check as well:
wc -l *35kb*txt > n_calls # get number of calls per chromosome
sed -i 's/^[ \t]*//' n_calls # fixes formatting of n_calls
grep '^0' n_calls > no_calls
sed -i 's/\./\t/g' no_calls
wc -l no_calls 
#0 no_calls
rm *n*calls

## Add chromosome names to each file, merge all chromosomes for each individual, and replace the individual's sample number with the individual's original file name (from 03_vcf_sample_order.list). 
for f in `seq 1 508`; do sed -e s/NUMBER/$f/g run.07_attach_chrom_names.sh > g.$f.sh; sbatch --mem=16000 g.$f.sh; rm g.$f.sh; done
for h in `seq 1 508`; do cat $h.35kb.d2.* > $h.35kb.d2.txt; done # can also use run.07 run.07_concatenate_chrom_files.sh
for h in `seq 1 508`; do tmp=`head -$h 04_vcf_sample_order.list | tail -1`; mv $h.35kb.d2.txt $tmp.35kb.d2.txt; done
## Add column with the individual ID
for f in `cat 04_vcf_sample_order.list`; do awk 'BEGIN {OFS="\t"} {print $0,FILENAME}' $f.35kb.d2.*txt > tmp.$f.txt; mv tmp.$f.txt $f.35kb.d2.txt; sed -i 's/.35kb.d2.txt//g' $f.35kb.d2.txt; done # can also use sbatch run.07b.add_indiv_ids.sh 

## Clean up step: get rid of genotype likelihoods, r.06_*, calls by chromosome, and any other intermediate files. 

## For each individual, [Tauras, not sure what you were planning on adding here?]
## For each ancestry informative SNP, use majority rule to assign the likely ancestry state using all SNPs within 35kb of that site.
## Note: majority rule requires that that call be at least 50% of calls in that window; no call will be made for a site without consensus. 
## Form ancestry tracts by combining continguous SNPs that have the same ancestry assignment
## Remove ancestry tracts of < 1kb, and output again. 
## Must module load R first before running run.08.mode_tracts.R script --> R/3.6.1-gcb03
for f in `cat 04_vcf_sample_order.list`; do cat run.08.mode_tracts.R | sed -e s/INDIV/$f/g > g.$f.sh; sbatch --mem=45000 g.$f.sh; done
## clean up the tmp2.INDIV.txt files, leaving just the majority rule files and the tracts. 

for f in `awk 'FNR>200 && FNR<=350' 04_vcf_sample_order.list`; do cat run.08.mode_tracts.R | sed -e s/INDIV/$f/g > g.$f.sh; sbatch --mem=45000 g.$f.sh; done

## Merge all calls into a single file
touch all.1kb.tracts.fullref.txt; for f in `ls *1kb.d2.*`; do sed '1d' $f >> all.1kb.tracts.fullref.txt; done

# for fig 1,
touch all.majrule.maskedSWref_and_unmaskedref.txt; for f in `ls *Maj*`; do sed '1d' $f >>  all.majrule.maskedSWref_and_unmaskedref.txt; done

for f in `cat 00_amboseli.list`; do mv $f*.d2*txt ancestry_usingfullref/amboseli_indiv/; done
for f in `cat 00_allref.list`; do mv $f*d2*txt ancestry_usingfullref/refpanel_indiv/; done

# In amboseli_indiv directory, merge majority rule call and tract files
touch amb.tracts.txt; for f in `ls *tracts.d2.*`; do sed '1d' $f >> amb.tracts.txt; done
touch amb.majrule.txt; for f in `ls *Maj*`; do sed '1d' $f >> amb.majrule.txt; done
