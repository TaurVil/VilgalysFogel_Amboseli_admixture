## _F4_ statistics

In addition to estimating admixture using a local ancestry approach (e.g., LCLAE), we used F4-ratio estimation (see Patterson et al. 2012 Genetics for details) as an orthogonal method to estimate the anubis ancestry proportions of Amboseli individuals. Following the nomenclature of Patterson et al. (2012), we estimated the ancestry proportion of members of population X (Amboseli) as a combination of two sources, represented by modern populations B (the SNPRC anubis founders; n=24) and (yellow baboons from Mikumi National Park in Tanzania; n=10). This method also requires two additional populations: population A, which has not contributed to population X but is a sister group to population B (hamadryas baboons and Guinea baboons resequenced by the Baboon Genome Sequencing Consortium as two alternatives; n=2 in each case), and an outgroup, O, to populations A, B, and C (the gelada monkey, also from the Baboon Genome Sequencing Consortium, and rhesus macaque; n=1 in each case).

Details can be found in Supplementary Methods Section 8.

XXX Once ASF gets SRR file names for new data, can update a few files below XXXX

```console

# Data needed:
# Fastq files for Amboseli high coverage individuals (recently admixed animals: AMB_301, AMB_317; historically admixed animals: AMB_310, AMB_311, AMB_312, AMB_313, AMB_314, AMB_316, AMB_318) are available in the NCBI Sequence Read Archive (SRA) under BioProject accession numbers PRJNAxxx and PRJNA295782 (see Table S1 for accession numbers for each sample)
# Fastq files for the SNPRC anubis founders are available in the NCBI Sequence Read Archive (SRA) under BioProject accession number PRJNA433868 (see Table S1 for accession numbers for each sample)
# Fastq files for yellow baboons from Mikumi National Park are available in the NCBI Sequence Read Archive (SRA) under BioProject accession number PRJNAxxx (see Table S1 for accession numbers for each sample)
# Fastq files for hamadryas, Guinea, and the gelada monkey from the Baboon Genome Sequencing Consortium available in the NCBI Sequence Read Archive (SRA) under BioProject accession numbers PRJNA20425, PRJNA54003, PRJNA251424, respectively (see Table S1 for accession numbers for each sample).

# Map these data to the rhesus macaque genome (MacaM; Zimin et al. 2014 Biology Direct) using bowtie2
for f in `cat XXXOnce you get SRR numbersXXX`; do sed -e s/NAME/$f/g 01a.bowtie2_map_macam.SRRs1.sh > $f.sh; sbatch --mem=32000 --out=mapping.$f.out $f.sh; done; rm $f.sh # for high coverage Amboseli data (excluding data from NCBI SRA BioProject PRJNA295782) and Mikumi data where paired-end reads are labeled as R1 and R2

cat SRR_list_BGDP SRR_list_SNPRC_amboseli_published SRR_list_SNPRC_anubis_founders >> SRR_tmp
for f in `cat SRR_tmp`; do sed -e s/NAME/$f/g 01a.bowtie2_map_macam.SRRs2.sh > $f.sh; sbatch --mem=32000 --out=mapping.$f.out $f.sh; done; rm $f.sh; rm SRR_tmp # for all other data where paired-end reads are labeled as _1 and _2

# Make directory to store mapped bams and move all these files to this directory
mkdir bams
mv map* bams/

# Merge multiple bams from the same individual into one bam per individual
# Must first sort before merging
for f in `cat SRRs_AMB_310`; do sed -e s/NAME/$f/g 02a.sort_for_merge.sh > $f.sh; sbatch --mem=2G $f.sh; done; rm $f.sh
for f in `cat SRRs_1X1126`; do sed -e s/NAME/$f/g 02a.sort_for_merge.sh > $f.sh; sbatch --mem=2G $f.sh; done; rm $f.sh
for f in `cat SRRs_1X1765`; do sed -e s/NAME/$f/g 02a.sort_for_merge.sh > $f.sh; sbatch --mem=2G $f.sh; done; rm $f.sh
for f in `cat SRR_list_BGDP`; do sed -e s/NAME/$f/g 02a.sort_for_merge.sh > $f.sh; sbatch --mem=2G $f.sh; done; rm $f.sh

# For individuals with two bams to merge (listed in fastqs_per_indiv2 where the first column is the individual's id which will be used as the new name of the merged file and the remaining columns are the bams to be merged), run:
for i in `seq 1 2`; do sed "${i}q;d" fastqs_per_indiv2 > tmp2; f=`awk '{print $1}' tmp2`; g=`awk '{print $2}' tmp2`; h=`awk '{print $3}' tmp2`; sed -e s/SAMPLE1/$f/g 02b.merge_2.sh > g.$f.sh; sed -e s/SAMPLE2/$g/g g.$f.sh > g.$f.$g.sh; sed -e s/SAMPLE3/$h/g  g.$f.$g.sh > 2.$f.$g.$h.sh; sbatch --mem=10G 2.$f.$g.$h.sh; done; rm 2.$f.$g.$h.sh

# For individuals with three bams to merge (listed in fastqs_per_indiv3 where the first column is the individual's id which will be used as the new name of the merged file and the remaining columns are the bams to be merged), run:
for i in `seq 1 3`; do sed "${i}q;d" fastqs_per_indiv3 > tmp2; f=`awk '{print $1}' tmp2`; g=`awk '{print $2}' tmp2`; h=`awk '{print $3}' tmp2`; k=`awk '{print $4}' tmp2`; sed -e s/SAMPLE1/$f/g 02b.merge_3.sh > g.$f.sh; sed -e s/SAMPLE2/$g/g g.$f.sh > g.$f.$g.sh;  sed -e s/SAMPLE3/$h/g  g.$f.$g.sh >  g.$f.$g.$h.sh;  sed -e s/SAMPLE4/$k/g g.$f.$g.$h.sh > 3.$f.$g.$h.$k.sh; sbatch --mem=10G 3.$f.$g.$h.$k.sh; done; rm 3.$f.$g.$h.$k.sh # note that we added the species name "Gelada" in front of the id # in fastqs_per_indiv3 to help distinguish individuals from different species and because all other ids from the Baboon Sequencing Genome Consortion have species as part of the id

# For individuals with four bams to merge (listed in fastqs_per_indiv4 where the first column is the individual's id which will be used as the new name of the merged file and the remaining columns are the bams to be merged), run:
for i in `seq 1 4`; do sed "${i}q;d" fastqs_per_indiv4 > tmp2; f=`awk '{print $1}' tmp2`; g=`awk '{print $2}' tmp2`; h=`awk '{print $3}' tmp2`; j=`awk '{print $4}' tmp2`; k=`awk '{print $5}' tmp2`; sed -e s/SAMPLE1/$f/g 02b.merge_4.sh > g.$f.sh; sed -e s/SAMPLE2/$g/g g.$f.sh > g.$f.$g.sh;  sed -e s/SAMPLE3/$h/g  g.$f.$g.sh >  g.$f.$g.$h.sh; sed -e s/SAMPLE4/$j/g g.$f.$g.$h.sh > g.$f.$g.$h.$j.sh;  sed -e s/SAMPLE5/$k/g g.$f.$g.$h.$j.sh > 4.$f.$g.$h.$j.$k.sh; sbatch --mem=10G  4.$f.$g.$h.$j.$k.sh; done; rm 4.$f.$g.$h.$j.$k.sh

# For individuals with five bams to merge (listed in fastqs_per_indiv5 where the first column is the individual's id which will be used as the new name of the merged file and the remaining columns are the bams to be merged), run:
for i in `seq 1 5`; do sed "${i}q;d" fastqs_per_indiv5 > tmp2; f=`awk '{print $1}' tmp2`; g=`awk '{print $2}' tmp2`; h=`awk '{print $3}' tmp2`; j=`awk '{print $4}' tmp2`; k=`awk '{print $5}' tmp2`; l=`awk '{print $6}' tmp2`; sed -e s/SAMPLE1/$f/g 02b.merge_5.sh > g.$f.sh; sed -e s/SAMPLE2/$g/g g.$f.sh > g.$f.$g.sh;  sed -e s/SAMPLE3/$h/g  g.$f.$g.sh >  g.$f.$g.$h.sh; sed -e s/SAMPLE4/$j/g g.$f.$g.$h.sh > g.$f.$g.$h.$j.sh;  sed -e s/SAMPLE5/$k/g g.$f.$g.$h.$j.sh > g.$f.$g.$h.$j.$k.sh; sed -e s/SAMPLE6/$l/g g.$f.$g.$h.$j.$k.sh > 5.$f.$g.$h.$j.$k.$l.sh; sbatch --mem=20G 5.$f.$g.$h.$j.$k.$l.sh;   done; rm 5.$f.$g.$h.$j.$k.$l.sh

mv map* bams/

# Make new list including all samples we will analyze
cat SRR_list* >> tmp # grab list of all SRR file names
cat SRRs* >> tmp2 # grab list of all SRR file names that were merged
awk 'NR==FNR{a[$0];next} !($0 in a)' tmp2 tmp >> tmp3 # remove all SRR file names that were merged per individual
awk '{print $1}' fastqs_per_indiv* >> tmp4
cat tmp3 tmp4 >> samples_final # add merged bam files per individual to sample list
rm tmp* # remove temporary files that are no longer needed

# Process samples before generating gVCF files
# Sort, add read groups, remove duplicate reads and reads with MAPQ less than 10
for f in `cat samples_final`; do sed -e s/SAMPLE/$f/g 02c.sort_RG_nodup_mapq10.sh > g.sh; sbatch --mem=30000 g.sh; done; rm g.sh 

# Make directory to store gVCFs 
mkdir gVCF

# Generate gVCF files using GATK HaplotypeCaller requiring a minimum base quality ≥ 20
# Run for each chromosome and each individual
for f in `cat MacaM_autosome_list`; do sed -e s/CHROM/$f/g 03.gvcf.sh > $f.sh; sbatch --array=1-50 --mem=15000 $f.sh; done; rm $f.sh
# Merge gVCF files across individuals using GATK CombineGVCF and then call genotypes using GATK GenotypeGVCFs. Also, filter for high quality variants following GATK’s recommended criteria for hard filtering for germline variants and retain biallelic SNPs that were typed within all individuals in the sample. In addition, remove indels, singleton and doubleton variants, and clusters of 3 or more variants that fall within a 10 bp window.
for f in `cat MacaM_autosome_list`; do sed -e s/CHROM/$f/g r04.merge_gvcfs_genotype_filter.sh > $f.sh; sbatch --mem=30000 $f.sh; done; rm $f.sh

# We would like to use the macaque as one possible outgroup so add the macaque genotype (homozygous reference) as an additional sample at all sites
for f in `cat MacaM_autosome_list`; do sed -e s/CHROMOSOME/$f/g 05.add_macaque_geno.sh >> $f.sh; sbatch --mem=20 $f.sh; done; rm $f.sh

# Convert each chromosome-specific vcf into EIGENSTRAT format by way of the plink format
# Note that we will convert chromosomes 1-19 (excluding chromosomes 02a and 02b) separately from chromosomes 02a and 02b which require some name adjustments for plink
# This step requires the plink (https://www.cog-genomics.org/plink2) and the covertf program from EIGENSOFT (https://github.com/DReichLab/EIG; https://github.com/DReichLab/EIG/tree/master/CONVERTF)
# Before running the script to convert vcf --> plink --> EIGENSTRAT, create a 'parfile' which is required for the the plink --> EIGENSTRAT conversion using convertf 
for f in `cat MacaM_autosome_list`; do sed -e s/CHROM/$f/g my.par.ped.eigenstrat >> my.par.ped.eigenstrat.$f; done # create chromosome-specific my.par.ped.eigenstrat files; note that macaque is specified as the sample id (polarize=macaque) so that reference/alternate alleles are always determined by macaque which will have a count of 2 (note that if this is not specified, we'll stil get the same f4 results because it doesn't matter which is allele is called reference/alternative)

# For chromosomes 1-19 (excluding chromosomes 02a and 02b), run vcf --> plink --> EIGENSTRAT conversion:
grep -v "chr02" MacaM_autosome_list >> MacaM_autosome_list_no_chr02s # get list of autosomes excluding chr02a and chr02b
for f in `cat MacaM_autosome_list_no_chr02s`; do sed -e s/CHROM/$f/g 06a.get_eigenstrat_format_excluding_chr02s.sh >> $f.sh; sbatch $f.sh; done; rm MacaM_autosome_list_no_chr02s; rm $f.sh # for all autosomes except chr02a and chr02b

# For chromosomes 02a and 02b, plink and EIGENSTRAT formats do not like weird chromosome names so rename chr02a and chr02b to chr20 and chr21 in the vcf, respectively (but can keep file names as chr02a and chr02b)
echo "chr02a chr20" >> chr02a_names
echo "chr02b chr21" >> chr02b_names

module load bcftools
module load tabix
bcftools annotate --rename-chrs chr02a_names final.baboon.to.macam.n49.chr02a.vcf.gz | bgzip > final.baboon.to.macam.n49.chr02a.rename.vcf.gz
bcftools annotate --rename-chrs chr02b_names final.baboon.to.macam.n49.chr02b.vcf.gz | bgzip > final.baboon.to.macam.n49.chr02b.rename.vcf.gz
tabix final.baboon.to.macam.n49.chr02a.rename.vcf.gz 
tabix final.baboon.to.macam.n49.chr02b.rename.vcf.gz 

rm chr02*names

# Now run vcf --> plink --> EIGENSTRAT conversion:
grep "chr02" MacaM_autosome_list >> MacaM_autosome_list_only_chr02s # get list of chr02a and chr02b
for f in `cat MacaM_autosome_list_only_chr02s`; do sed -e s/CHROM/$f/g 06a.get_eigenstrat_format_chr02s.sh >> $f.sh; sbatch $f.sh; done; rm MacaM_autosome_list_only_chr02s; rm $f.sh # for renamed chr02a and chr02b files

# The eig.n49.ind chromosome-specific files need to be updated with population info for each sample (see eig.indiv) and it's easier to remove the ones generated by convertf and copy eig.indiv so that each chromosome has the correct eig.ind file
rm eig.n48*ind # remove chromosome-specific files generated by convertf
for f in `cat MacaM_autosome_list`; do cp eig.indiv eig.n49.$f.ind; done # create new eig.ind files

# Now we have the requisite EIGENSTRAT files (ind, snp, and geno files for each chromosome)!
# We'll input the EIGENSTRAT formatted data into admixr (version 0.9.1; Petr et al. 2019 Bioinformatics; https://github.com/bodkan/admixr), an R package that builds on the ADMIXTOOLS software suite (Patterson et al. 2012 Genetics; https://github.com/DReichLab/AdmixTools). There are several helpful tutorials for using admixr (see here for example: https://bodkan.net/admixr/articles/tutorial.html).
# admixr is available on CRAN and be installed like a typical R package using:
install.packages("admixr")
# However, admixr requires that ADMIXTOOLS is also downloaded. After installing ADMIXTOOLS (https://github.com/DReichLab/AdmixTools/blob/master/README.INSTALL), you'll need to load a few more packages for ADMIXTOOLS to work:
module load gcc
module load gsl
module load OpenBLAS
# You'll also need to change your path so R can find where ADMIXTOOLS is installed
export PATH=$PATH:AdmixTools-master/bin
# Check that it's in the PATH
echo $PATH # should be the last part of the path

# Estimate the f4 ratio per chromosome using the f4ratio function in admixr
for f in `cat MacaM_autosome_list`; do sed -e s/CHROM/$f/g 07.f4ratio_admixr.R >> $f.sh; sbatch $f.sh; done; rm $f.sh

# To obtain a genome-wide f4 ratio estimate of anubis ancestry proportions (α), averaged for the different phylogenetic configurations per chromosome and then averaged across chromosomes, weighted by chromosome length, run:
8.f4ratio_result.R

```
