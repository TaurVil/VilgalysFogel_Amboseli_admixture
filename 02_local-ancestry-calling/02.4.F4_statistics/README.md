## _F4_ statistics

In addition to estimating admixture using a local ancestry approach (e.g., LCLAE), we used _F4_-ratio estimation (see Patterson et al. 2012 _Genetics_ for details) as an orthogonal method to estimate the anubis ancestry proportions of Amboseli individuals. Following the nomenclature of Patterson et al. (2012), we estimated the ancestry proportion of members of population _X_ (Amboseli) as a combination of two sources, represented by modern populations _B_ (the SNPRC anubis founders; n=24) and _C_ (yellow baboons from Mikumi National Park in Tanzania; n=10). This method also requires two additional populations: population _A_, which has not contributed to population _X_ but is a sister group to population _B_ (hamadryas baboons and Guinea baboons resequenced by the Baboon Genome Sequencing Consortium as two alternatives; n=2 in each case), and an outgroup, _O_, to populations _A_, _B_, and _C_ (the gelada monkey, also from the Baboon Genome Sequencing Consortium, and rhesus macaque; n=1 in each case).

Details can be found in Supplementary Methods Section 8.

#### Mapping samples to the _MacaM_ genome

Obtain fastq files for high coverage individuals from Amboseli and unadmixed populations, then map these data to the rhesus macaque genome (_MacaM_; Zimin et al. 2014 Biology Direct) using bowtie2. Some of these files are split between multiple SRA ascensions, so they will need to be merged together into a single file after mapping. 

**Data needed:**
* Fastq files for Amboseli high coverage individuals (recently admixed animals: AMB_301, AMB_317; historically admixed animals: AMB_310, AMB_311, AMB_312, AMB_313, AMB_314, AMB_316, AMB_318) are available in the NCBI Sequence Read Archive (SRA) under BioProject accession PRJNA295782 or will be available at the time of publication under BioProject accession PRJNA755322 (see Table S1 for accession numbers for each sample)
* Fastq files for the SNPRC anubis founders are available in the NCBI Sequence Read Archive (SRA) under BioProject accession number PRJNA433868 (see Table S1 for accession numbers for each sample)
* Fastq files for yellow baboons from Mikumi National Park are available in the NCBI Sequence Read Archive (SRA) under BioProject accession PRJNA308870 or will be available at the time of publication under BioProject accession PRJNA755322 (see Table S1 for accession numbers for each sample)
* Fastq files for hamadryas, Guinea, and the gelada monkey from the Baboon Genome Sequencing Consortium available in the NCBI Sequence Read Archive (SRA) under BioProject accession numbers PRJNA20425, PRJNA54003, PRJNA251424, respectively (see Table S1 for accession numbers for each sample).

```console
# Map to MacaM using bowtie2
# SRR accession numbers will be made available at the time of publication
# for f in `cat SRR_newly_published`; do sed -e s/NAME/$f/g run.01a.bowtie2_map_macam.SRRs1.sh >> g.$f.sh; sbatch --mem=32000 --out=mapping.$f.out g.$f.sh; done # for high coverage Amboseli data (excluding data from NCBI SRA BioProject PRJNA295782) and Mikumi data where paired-end reads are labeled as R1 and R2
rm g*sh

# Map previously published data available on NCBI
cat SRR_list_BGDP SRR_list_SNPRC_amboseli_published SRR_list_SNPRC_anubis_founders >> SRR_tmp
for f in `cat SRR_tmp`; do sed -e s/NAME/$f/g run.01b.bowtie2_map_macam.SRRs2.sh >> g.$f.sh; sbatch --mem=32000 --out=mapping.$f.out g.$f.sh; done # for all other data where paired-end reads are labeled as _1 and _2
rm g*sh

# Make directory to store mapped bams and move all bam files to this directory
mkdir bams; mv map* bams/; rm SRR_tmp

# Merge multiple bams from the same individual into one bam per individual
# Must first sort before merging
for f in `cat SRRs_AMB_310`; do sed -e s/NAME/$f/g run.02a.sort_for_merge.sh >> g.$f.sh; sbatch --mem=2G g.$f.sh; done
for f in `cat SRRs_1X1126`; do sed -e s/NAME/$f/g run.02a.sort_for_merge.sh >> g.$f.sh; sbatch --mem=2G g.$f.sh; done
for f in `cat SRRs_1X1765`; do sed -e s/NAME/$f/g run.02a.sort_for_merge.sh >> g.$f.sh; sbatch --mem=2G g.$f.sh; done
for f in `cat SRR_list_BGDP`; do sed -e s/NAME/$f/g run.02a.sort_for_merge.sh >> g.$f.sh; sbatch --mem=2G g.$f.sh; done; rm g*sh

# For individuals with two bams to merge (listed in fastqs_per_indiv2 where the first column is the individual's id which will be used as the new name of the merged file and the remaining columns are the bams to be merged), run:
for i in `seq 1 2`; do sed "${i}q;d" fastqs_per_indiv2 > tmp2; f=`awk '{print $1}' tmp2`; g=`awk '{print $2}' tmp2`; h=`awk '{print $3}' tmp2`; sed -e s/SAMPLE1/$f/g run.02b.merge_2.sh > g.$f.sh; sed -e s/SAMPLE2/$g/g g.$f.sh > g.$f.$g.sh; sed -e s/SAMPLE3/$h/g  g.$f.$g.sh > 2.$f.$g.$h.sh; sbatch --mem=10G 2.$f.$g.$h.sh; done
rm 2.*.sh

# For individuals with three bams to merge (listed in fastqs_per_indiv3 where the first column is the individual's id which will be used as the new name of the merged file and the remaining columns are the bams to be merged), run:
for i in `seq 1 3`; do sed "${i}q;d" fastqs_per_indiv3 > tmp2; f=`awk '{print $1}' tmp2`; g=`awk '{print $2}' tmp2`; h=`awk '{print $3}' tmp2`; k=`awk '{print $4}' tmp2`; sed -e s/SAMPLE1/$f/g run.02b.merge_3.sh > g.$f.sh; sed -e s/SAMPLE2/$g/g g.$f.sh > g.$f.$g.sh;  sed -e s/SAMPLE3/$h/g  g.$f.$g.sh >  g.$f.$g.$h.sh;  sed -e s/SAMPLE4/$k/g g.$f.$g.$h.sh > 3.$f.$g.$h.$k.sh; sbatch --mem=10G 3.$f.$g.$h.$k.sh; done # note that we added the species name "Gelada" in front of the id # in fastqs_per_indiv3 to help distinguish individuals from different species and because all other ids from the Baboon Sequencing Genome Consortion have species as part of their ids
rm 3.*.sh

# For individuals with four bams to merge (listed in fastqs_per_indiv4 where the first column is the individual's id which will be used as the new name of the merged file and the remaining columns are the bams to be merged), run:
for i in `seq 1 4`; do sed "${i}q;d" fastqs_per_indiv4 > tmp2; f=`awk '{print $1}' tmp2`; g=`awk '{print $2}' tmp2`; h=`awk '{print $3}' tmp2`; j=`awk '{print $4}' tmp2`; k=`awk '{print $5}' tmp2`; sed -e s/SAMPLE1/$f/g run.02b.merge_4.sh > g.$f.sh; sed -e s/SAMPLE2/$g/g g.$f.sh > g.$f.$g.sh;  sed -e s/SAMPLE3/$h/g  g.$f.$g.sh >  g.$f.$g.$h.sh; sed -e s/SAMPLE4/$j/g g.$f.$g.$h.sh > g.$f.$g.$h.$j.sh;  sed -e s/SAMPLE5/$k/g g.$f.$g.$h.$j.sh > 4.$f.$g.$h.$j.$k.sh; sbatch --mem=10G  4.$f.$g.$h.$j.$k.sh; done
rm 4.*.sh

# For individuals with five bams to merge (listed in fastqs_per_indiv5 where the first column is the individual's id which will be used as the new name of the merged file and the remaining columns are the bams to be merged), run:
for i in `seq 1 5`; do sed "${i}q;d" fastqs_per_indiv5 > tmp2; f=`awk '{print $1}' tmp2`; g=`awk '{print $2}' tmp2`; h=`awk '{print $3}' tmp2`; j=`awk '{print $4}' tmp2`; k=`awk '{print $5}' tmp2`; l=`awk '{print $6}' tmp2`; sed -e s/SAMPLE1/$f/g run.02b.merge_5.sh > g.$f.sh; sed -e s/SAMPLE2/$g/g g.$f.sh > g.$f.$g.sh;  sed -e s/SAMPLE3/$h/g  g.$f.$g.sh >  g.$f.$g.$h.sh; sed -e s/SAMPLE4/$j/g g.$f.$g.$h.sh > g.$f.$g.$h.$j.sh;  sed -e s/SAMPLE5/$k/g g.$f.$g.$h.$j.sh > g.$f.$g.$h.$j.$k.sh; sed -e s/SAMPLE6/$l/g g.$f.$g.$h.$j.$k.sh > 5.$f.$g.$h.$j.$k.$l.sh; sbatch --mem=20G 5.$f.$g.$h.$j.$k.$l.sh; done
rm 5.*.sh

# Move merged bams to `bams` directory
mv map* bams/
```

For each individual, call genotypes using GATK. 

```console
# Make new list including all of the samples we will analyze
cat SRR_list* >> tmp # grab list of all SRR file names
cat SRRs* >> tmp2 # grab list of all SRR file names that were merged
awk 'NR==FNR{a[$0];next} !($0 in a)' tmp2 tmp >> tmp3 # remove all SRR file names that were merged because there were multiple files per individual
awk '{print $1}' fastqs_per_indiv* >> tmp4
cat tmp3 tmp4 >> samples_final # add merged bam files per individual to sample list
rm tmp* # remove temporary files that are no longer needed

# Process samples before generating gVCF files: sort, add read groups, remove duplicate reads, and reads with MAPQ less than 10
for f in `cat samples_final`; do sed -e s/SAMPLE/$f/g run.02c.sort_RG_nodup_mapq10.sh >> g.$f.sh; sbatch --mem=30000 g.$f.sh; done
rm g*sh

# Make directory to store gVCFs 
mkdir gVCF

# Generate gVCF files using GATK HaplotypeCaller requiring a minimum base quality ≥ 20
# Run for each chromosome and each individual
for f in `cat MacaM_autosome_list`; do sed -e s/CHROM/$f/g run.03.gvcf.sh >> g.$f.sh; sbatch --array=1-50 --mem=15000 g.$f.sh; done
rm g*sh
# Merge gVCF files across individuals using GATK CombineGVCF and then call genotypes using GATK GenotypeGVCFs. Also, filter for high quality variants following GATK’s recommended criteria for hard filtering for germline variants and retain biallelic SNPs that were typed within all individuals in the sample. In addition, remove indels, singleton and doubleton variants, and clusters of 3 or more variants that fall within a 10 bp window.
for f in `cat MacaM_autosome_list`; do sed -e s/CHROM/$f/g run.04.merge_gvcfs_genotype_filter.sh > g.$f.sh; sbatch --mem=30000 g.$f.sh; done
rm g*sh
# note: we performed joint genotyping with other samples that were subsequently removed before downstream analysis detailed in this directory. We do not anticipate that the inclusion of other samples during the joint genotyping step affects our results given our subsequent filtering and processing before analysis.

# to use the macaque as a possible outgroup, add the macaque genotype (homozygous reference) as an additional sample at all sites
for f in `cat MacaM_autosome_list`; do sed -e s/CHROMOSOME/$f/g run.05.add_macaque_geno.sh >> g.$f.sh; sbatch --mem=20 g.$f.sh; done
rm g*sh
```

Having completed joint genotype calling, we next need to convert each vcf into EIGENSTRAT format. We will then calculate _F4_ statistics using _admixr_ (Petr et al. 2019 _Bioinformatics_; https://github.com/bodkan/admixr)

```console
# Convert each chromosome-specific vcf into EIGENSTRAT format by way of the plink format
# Note that we will convert chromosomes 1-19 (excluding chromosomes 02a and 02b) separately from chromosomes 02a and 02b which require some name adjustments for plink/EIGENSTRAT
# This step requires the plink (https://www.cog-genomics.org/plink2) and the covertf program from EIGENSOFT (https://github.com/DReichLab/EIG; https://github.com/DReichLab/EIG/tree/master/CONVERTF)
# Before running the script to convert vcf --> plink --> EIGENSTRAT, create a 'parfile' which is required for the the plink --> EIGENSTRAT conversion using convertf 
for f in `cat MacaM_autosome_list`; do sed -e s/CHROM/$f/g my.par.ped.eigenstrat >> my.par.ped.eigenstrat.$f; done # create chromosome-specific my.par.ped.eigenstrat files; note that macaque is specified as the sample id (polarize=macaque) so that reference/alternate alleles are always determined by macaque which will have a count of 2 (note that if this is not specified, we'll stil get the same f4 results because it doesn't matter which allele is called reference/alternative)

# For chromosomes 1-19 (excluding chromosomes 02a and 02b), run vcf --> plink --> EIGENSTRAT conversion:
grep -v "chr02" MacaM_autosome_list >> MacaM_autosome_list_no_chr02s # get list of autosomes excluding chr02a and chr02b
for f in `cat MacaM_autosome_list_no_chr02s`; do sed -e s/CHROM/$f/g run.06a.get_eigenstrat_format_excluding_chr02s.sh >> g.$f.sh; sbatch g.$f.sh; done # for all autosomes except chr02a and chr02b
rm MacaM_autosome_list_no_chr02s; rm g*sh 

# For chromosomes 02a and 02b, plink and EIGENSTRAT formats do not like weird chromosome names so rename chr02a and chr02b to chr20 and chr21 in the vcf, respectively (but can keep file names as chr02a and chr02b)
echo "chr02a chr20" >> chr02a_names; echo "chr02b chr21" >> chr02b_names

module load bcftools; module load tabix
bcftools annotate --rename-chrs chr02a_names final.baboon.to.macam.n49.chr02a.vcf.gz | bgzip > final.baboon.to.macam.n49.chr02a.rename.vcf.gz
bcftools annotate --rename-chrs chr02b_names final.baboon.to.macam.n49.chr02b.vcf.gz | bgzip > final.baboon.to.macam.n49.chr02b.rename.vcf.gz
tabix final.baboon.to.macam.n49.chr02a.rename.vcf.gz; tabix final.baboon.to.macam.n49.chr02b.rename.vcf.gz 
rm chr02*names

# Now run vcf --> plink --> EIGENSTRAT conversion for chromosomes 02a and 02b:
grep "chr02" MacaM_autosome_list >> MacaM_autosome_list_only_chr02s # get list of chr02a and chr02b
for f in `cat MacaM_autosome_list_only_chr02s`; do sed -e s/CHROM/$f/g run.06b.get_eigenstrat_format_chr02s.sh >> g.$f.sh; sbatch g.$f.sh; done # for renamed chr02a and chr02b files
rm MacaM_autosome_list_only_chr02s; rm g*sh 

# Remove intermediate files no longer needed
rm updated*; rm *map; rm *ped

# The eig.n49.ind chromosome-specific files generated by the `run.06` scripts need to be updated with population info for each sample (see DATA/eig.indiv) and it's easier to remove the ones generated by convertf and copy eig.indiv so that each chromosome has the correct eig.ind file
rm eig.n49*ind # remove chromosome-specific "ind" files generated by convertf
for f in `cat MacaM_autosome_list`; do cp eig.indiv eig.n49.$f.ind; done # create new eig.ind files for each chromosome

# Now we have the requisite EIGENSTRAT files (ind, snp, and geno files for each chromosome)!
# We'll input the EIGENSTRAT formatted data into admixr (version 0.9.1; Petr et al. 2019 Bioinformatics; https://github.com/bodkan/admixr), an R package that builds on the ADMIXTOOLS software suite (Patterson et al. 2012 Genetics; https://github.com/DReichLab/AdmixTools). There are several helpful tutorials for using admixr (see here for an example: https://bodkan.net/admixr/articles/tutorial.html).
# admixr is available on CRAN and be installed like a typical R package using:
install.packages("admixr")
# However, admixr requires that ADMIXTOOLS is also downloaded. After installing ADMIXTOOLS (https://github.com/DReichLab/AdmixTools/blob/master/README.INSTALL), you'll need to load a few more packages for ADMIXTOOLS to work:
module load gcc
module load gsl
module load OpenBLAS
# You'll also need to change your path so R can find where ADMIXTOOLS is installed
export PATH=$PATH:AdmixTools-master/bin # location of the AdmixTools-master directory
# Check that it's in the PATH
echo $PATH # should be the last part of the path

# Estimate the f4 ratio per chromosome using the f4ratio function in admixr
for f in `cat MacaM_autosome_list`; do sed -e s/CHROM/$f/g run.07.f4ratio_admixr.R >> g.$f.sh; sbatch g.$f.sh; done
rm g*sh

# To obtain a genome-wide f4 ratio estimate of anubis ancestry proportions (α), get the average of the different phylogenetic configurations per chromosome and then average these estimates across chromosomes, weighted by chromosome length. We provide the chromosome-specific Rdata files generated by the run.07.f4ratio_admixr.R script in DATA (titled f4_ratio_CHROMOSOME.Rd) in order to run the following script. Run the R script:
run.08.f4ratio_results.R
```
