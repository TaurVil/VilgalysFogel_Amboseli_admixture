## F4 statistics

In addition to estimating admixture using a local ancestry approach (e.g., LCLAE), we used F4-ratio estimation (see Patterson et al. 2012 Genetics for details) as an orthogonal method to estimate the anubis ancestry proportions of Amboseli individuals. Following the nomenclature of Patterson et al. (2012), we estimated the ancestry proportion of members of population X (Amboseli) as a combination of two sources, represented by modern populations B (the SNPRC anubis founders; n=24) and (yellow baboons from Mikumi National Park in Tanzania; n=10). This method also requires two additional populations: population A, which has not contributed to population X but is a sister group to population B (hamadryas baboons and Guinea baboons resequenced by the Baboon Genome Sequencing Consortium as two alternatives; n=2 in each case), and an outgroup, O, to populations A, B, and C (the gelada monkey, also from the Baboon Genome Sequencing Consortium, and rhesus macaque; n=1 in each case). 

Details can be found in Supplementary Methods Section 8.

```console

# Data needed:
# Fastq files for Amboseli high coverage individuals (recently admixed animals: AMB_301, AMB_317; historically admixed animals: AMB_310, AMB_311, AMB_312, AMB_313, AMB_314, AMB_316, AMB_318) are available in the NCBI Sequence Read Archive (SRA) under BioProject accession numbers PRJNAxxx and PRJNA295782 (see Table S1 for accession numbers for each sample)
# Fastq files for the SNPRC anubis founders are available in the NCBI Sequence Read Archive (SRA) under BioProject accession number PRJNA433868 (see Table S1 for accession numbers for each sample)
# Fastq files for yellow baboons from Mikumi National Park are available in the NCBI Sequence Read Archive (SRA) under BioProject accession number PRJNAxxx (see Table S1 for accession numbers for each sample)
# Fastq files for hamadryas, Guinea, and the gelada monkey from the Baboon Genome Sequencing Consortium available in the NCBI Sequence Read Archive (SRA) under BioProject accession numbers PRJNA20425, PRJNA54003, PRJNA251424, respectively (see Table S1 for accession numbers for each sample).

# Map these data to the rhesus macaque genome (MacaM; Zimin et al. 2014 Biology Direct) using bowtie2
for f in `cat XXXOnce you get SRR numbersXXX`; do sed -e s/NAME/$f/g 01a.bowtie2_map_macam.SRRs1.sh > $f.sh; sbatch --mem=32000 --out=mapping.$f.out $f.sh; rm $f.sh; done # for high coverage Amboseli data (excluding data from NCBI SRA BioProject PRJNA295782) and Mikumi data where paired-end reads are labeled as R1 and R2

cat SRR_list_BGDP SRR_list_SNPRC_amboseli_published SRR_list_SNPRC_anubis_founders >> SRR_tmp
for f in `cat SRR_tmp`; do sed -e s/NAME/$f/g 01a.bowtie2_map_macam.SRRs2.sh > $f.sh; sbatch --mem=32000 --out=mapping.$f.out $f.sh; rm $f.sh; rm SRR_tmp; done # for all other data where paired-end reads are labeled as _1 and _2

# Make directory to store mapped bams and move all these files to this directory
mkdir bams
mv map* bams/

# Merge multiple bams from the same individual into one bam per individual
# Must first sort before merging
for f in `cat LIT_list`; do sed -e s/NAME/$f/g 02a.sort_for_merge.sh > $f.sh; sbatch --mem=2G $f.sh; rm $f.sh; done 
for f in `cat 1X1765_list_list`; do sed -e s/NAME/$f/g 02a.sort_for_merge.sh > $f.sh; sbatch --mem=2G $f.sh; rm $f.sh; done
for f in `cat 1X1765_list`; do sed -e s/NAME/$f/g 02a.sort_for_merge.sh > $f.sh; sbatch --mem=2G $f.sh; rm $f.sh; done 
for f in `cat 00_BGDP`; do sed -e s/NAME/$f/g 02a.sort_for_merge.sh > $f.sh; sbatch --mem=2G $f.sh; rm $f.sh; done

# rename samples with only one FASTQ
for i in `seq 1 12`; do sed "${i}q;d" fastqs_per_indiv1 > tmp2; f=`awk '{print $1}' tmp2`; g=`awk '{print $2}' tmp2`; sed -e s/SAMPLE1/$f/g run.02b.rename_1.sh | sed -e s/SAMPLE/$g/g > 1.$f.$g.sh; sbatch 1.$f.$g.sh; done

# merge 3 samples from the same individual
for i in `seq 1`; do sed "${i}q;d" fastqs_per_indiv3 > tmp2; f=`awk '{print $1}' tmp2`; g=`awk '{print $2}' tmp2`; h=`awk '{print $3}' tmp2`; k=`awk '{print $4}' tmp2`; sed -e s/SAMPLE1/$f/g run.02b.merge_3.sh > g.$f.sh; sed -e s/SAMPLE2/$g/g g.$f.sh > g.$f.$g.sh;  sed -e s/SAMPLE3/$h/g  g.$f.$g.sh >  g.$f.$g.$h.sh;  sed -e s/SAMPLE4/$k/g g.$f.$g.$h.sh > 3.$f.$g.$h.$k.sh; sbatch --mem=10G 3.$f.$g.$h.$k.sh; done
# merge 4 samples from the same individual
for i in `seq 1 3`; do sed "${i}q;d" fastqs_per_indiv4 > tmp2; f=`awk '{print $1}' tmp2`; g=`awk '{print $2}' tmp2`; h=`awk '{print $3}' tmp2`; j=`awk '{print $4}' tmp2`; k=`awk '{print $5}' tmp2`; sed -e s/SAMPLE1/$f/g run.02b.merge_4.sh > g.$f.sh; sed -e s/SAMPLE2/$g/g g.$f.sh > g.$f.$g.sh;  sed -e s/SAMPLE3/$h/g  g.$f.$g.sh >  g.$f.$g.$h.sh; sed -e s/SAMPLE4/$j/g g.$f.$g.$h.sh > g.$f.$g.$h.$j.sh;  sed -e s/SAMPLE5/$k/g g.$f.$g.$h.$j.sh > 4.$f.$g.$h.$j.$k.sh; sbatch --mem=10G  4.$f.$g.$h.$j.$k.sh; done
# merge 5 samples from the same individual
for i in `seq 1 5`; do sed "${i}q;d" fastqs_per_indiv5 > tmp2; f=`awk '{print $1}' tmp2`; g=`awk '{print $2}' tmp2`; h=`awk '{print $3}' tmp2`; j=`awk '{print $4}' tmp2`; k=`awk '{print $5}' tmp2`; l=`awk '{print $6}' tmp2`; sed -e s/SAMPLE1/$f/g run.02b.merge_5.sh > g.$f.sh; sed -e s/SAMPLE2/$g/g g.$f.sh > g.$f.$g.sh;  sed -e s/SAMPLE3/$h/g  g.$f.$g.sh >  g.$f.$g.$h.sh; sed -e s/SAMPLE4/$j/g g.$f.$g.$h.sh > g.$f.$g.$h.$j.sh;  sed -e s/SAMPLE5/$k/g g.$f.$g.$h.$j.sh > g.$f.$g.$h.$j.$k.sh; sed -e s/SAMPLE6/$l/g g.$f.$g.$h.$j.$k.sh > 5.$f.$g.$h.$j.$k.$l.sh; sbatch --mem=20G  5.$f.$g.$h.$j.$k.$l.sh; done
# merge 8 samples from the same individual
for i in `seq 1`; do sed "${i}q;d" fastqs_per_indiv8 > tmp2; f=`awk '{print $1}' tmp2`; g=`awk '{print $2}' tmp2`; h=`awk '{print $3}' tmp2`; j=`awk '{print $4}' tmp2`; k=`awk '{print $5}' tmp2`; l=`awk '{print $6}' tmp2`; m=`awk '{print $7}' tmp2`; n=`awk '{print $8}' tmp2`; o=`awk '{print $9}' tmp2`; sed -e s/SAMPLE1/$f/g run.02b.merge_5.sh > g.$f.sh; sed -e s/SAMPLE2/$g/g g.$f.sh > g.$f.$g.sh;  sed -e s/SAMPLE3/$h/g  g.$f.$g.sh >  g.$f.$g.$h.sh; sed -e s/SAMPLE4/$j/g g.$f.$g.$h.sh > g.$f.$g.$h.$j.sh;  sed -e s/SAMPLE5/$k/g g.$f.$g.$h.$j.sh > g.$f.$g.$h.$j.$k.sh; sed -e s/SAMPLE6/$l/g g.$f.$g.$h.$j.$k.sh > g.$f.$g.$h.$j.$k.$l.sh; sed -e s/SAMPLE7/$m/g g.$f.$g.$h.$j.$k.$l.sh > g.$f.$g.$h.$j.$k.$l.$m.sh; sed -e s/SAMPLE8/$n/g g.$f.$g.$h.$j.$k.$l.$m.sh > g.$f.$g.$h.$j.$k.$l.$m.$n.sh; sed -e s/SAMPLE9/$o/g g.$f.$g.$h.$j.$k.$l.$m.$n.sh > 8.$f.$g.$h.$j.$k.$l.$m.$n.$o.sh; sbatch --mem=20G 8.$f.$g.$h.$j.$k.$l.$m.$n.$o.sh; done

sbatch --mem=1G run.02b.merge_LIT.sh
sbatch --mem=1G run.02b.merge_SW_1X1765.sh 
sbatch --mem=1G run.02b.merge_SW_1X1126.sh

# get individuals to merge for BGDP
# in R
library(dplyr)
tmp <- read.table("SraRunTable_BGDP_nonyelanu", header=T, sep=",")
tmp$Sample.Name2 <- gsub("-", ".", tmp$Sample.Name)
# add Gelada species name to sample (all other sample names have their species in their sample name)
which(tmp$Sample.Name2=="38168")
tmp[54:56,]$Sample.Name2 <- "GELADA.38168"
tmp2 <- tmp %>% group_by(Sample.Name2) %>% tally()
# do for samples with 1, 3, 4, 5, and 8 fastqs per individual separately (showing example for 8)
poop <- tmp2[tmp2$n==8,]
poop2 <- tmp[tmp$Sample.Name2 %in% poop$Sample.Name2,]
poop3 <- poop2[c(37,1)]
poop4 <- poop3 %>% group_by(Sample.Name2) %>% summarize(value=paste(Run, collapse="\t"))
write.table(poop4, "fastqs_per_indiv8", row.names=F, col.names=F, sep="\t", quote=F)



cat 00* >> 00_all_names
cat *list >> multiple_bams_list 
awk 'NR==FNR{a[$0];next} !($0 in a)' multiple_bams_list 00_all_names >> 00_all_names2
# add
#1X1126
#1X1765
#LIT
# to the end of the 00_all_names2

mv 1X1126.bam mapped.1X1126.bam
mv 1X1765.bam mapped.1X1765.bam
mv LIT.bam mapped.LIT.bam

for f in `head 00_all_names2 | tail -n +8`; do sed -e s/SAMPLE/$f/g run.02c.sort_RG_nodup_mapq10.sh > g.sh; sbatch --mem=30000 g.sh; done; rm g.sh 
for f in `cat merged_bam_list`; do sed -e s/SAMPLE/$f/g run.02c.sort_RG_nodup_mapq10.sh > g.sh; sbatch --mem=30000 g.sh; done; rm g.sh # bgdp samples



mkdir gVCF

for f in `cat MacaM_autosome_list`; do sed -e s/CHROM/$f/g run.03.gvcf.sh > $f.sh; sbatch --array=1-94%16 --mem=15000 $f.sh; done
for f in `cat MacaM_autosome_list`; do sed -e s/CHROM/$f/g run.04a.merge_gvcfs.sh > $f.sh; sbatch --mem=30000 $f.sh; done

sbatch --mem=12000 run.04a.merge_gvcfs.sh


# or run for f in `cat MacaM_autosome_list`; do sed -e s/CHROMOSOME/$f/g run.05.add_macaque_geno.sh >> $f.sh; sbatch --mem=20 $f.sh; done
rm chr*sh
# we now need to add the macaque genotypes - let's copy our vcf (named it tmp.vcf.gz) and try this out
module load bcftools
# remove vcf header
bcftools view --no-header biallelic.males.females.chrX.baboon.to.macam.recode.vcf >> noheader.biallelic.males.females.chrX.baboon.to.macam.recode.vcf
# add a new sample for all sites where the genotype is homozygous reference (I just grabbed the other values from one of the other reference genotypes - may want to change this later)
sed 's/$/\t0\/0:11,0:11:33:0,33,430/g' noheader.biallelic.males.females.chrX.baboon.to.macam.recode.vcf >> tmp.vcf

# add macaque as the name of our newly added sample
sed '/^#CHROM/ {s/$/\tmacaque/}' biallelic.males.females.chrX.baboon.to.macam.recode.vcf >> tmp2.vcf
# bgzip and index our new vcf
module load tabix
bgzip tmp2.vcf
tabix tmp2.vcf.gz
# only output the vcf header
bcftools view -h tmp2.vcf.gz >> header
###bcftools annotate tmp.vcf -h header
cat header tmp.vcf >> final.chrX.baboon.to.macam.recode.vcf

bgzip final.chrX.baboon.to.macam.recode.vcf 
tabix final.chrX.baboon.to.macam.recode.vcf.gz

# we now need to convert our vcf into eigenstrat format by way of the plink format
# plink/eigenstrat is does not like weird chromosome names so rename chr02a and chr02b to chr20 and chr21 respectively
echo "chr02a chr20" >> chr02a_names
echo "chr02b chr21" >> chr02b_names
module load bcftools
module load tabix
bcftools annotate --rename-chrs chr02a_names final.baboon.to.macam.n51.chr02a.vcf.gz | bgzip > final.baboon.to.macam.n51.chr02a.rename.vcf.gz
bcftools annotate --rename-chrs chr02b_names final.baboon.to.macam.n51.chr02b.vcf.gz | bgzip > final.baboon.to.macam.n51.chr02b.rename.vcf.gz
tabix final.baboon.to.macam.n51.chr02a.rename.vcf.gz 
tabix final.baboon.to.macam.n51.chr02b.rename.vcf.gz 

# first convert to plink, then add population info to plink ped file, then convert to eigenstrat
for f in `cat MacaM_autosome_list`; do sed -e s/CHROM/$f/g my.par.ped.eigenstrat >> my.par.ped.eigenstrat.$f; done # create chromosome specific my.par.ped.eigenstrat files
# note that macaque is specified as the sample id so that reference/alternate alleles are always determined by macaque which will have a count of 2 (note that if you don't specify this, you'll stil get the same f4 results because it doesn't matter which is allele is called reference/alternative)
for f in `cat MacaM_autosome_list_nochr2 `; do sed -e s/CHROM/$f/g run.05.get_eigenstrat_format.sh >> $f.sh; sbatch $f.sh; done # for all autosomes except chr02a and chr02b

grep "chr02" MacaM_autosome_list >> MacaM_autosome_list_chr02s
for f in `cat MacaM_autosome_list_chr02s`; do sed -e s/CHROM/$f/g run.05.get_eigenstrat_format_chr02s.sh >> $f.sh; sbatch $f.sh; done # for renamed chr02a and chr02b

# chromosome conversion between Rogers et al. 2006 Genomics linkage (see Table 1 linkage map) map/genome assembly and Zimin et al. 2014 Biol Direct (MacaM)
https://www.ncbi.nlm.nih.gov/genome/guide/rhesus_macaque/rhesuschrtable.html



# not used for for now --a2-allele --real-ref-alleles

# we now have the requisite eigenstrat files (ind, snp, and geno files)
module load gcc
module load gsl
module load OpenBLAS
module load R
export PATH=$PATH:/data/tunglab/asf40/wgs_data/MedGenome_ftp/f4stats/AdmixTools-master/bin
# check that it's in the PATH
echo $PATH
#/nfs/software/helmod/apps/Core/gsl/1.16-fasrc02/bin:/nfs/software/helmod/apps/Core/netcdf/4.3.2-fasrc03/bin:/nfs/software/helmod/apps/Core/hdf5/1.8.12-fasrc04/bin:/nfs/software/helmod/apps/Core/R_core/3.6.1-gcb03/bin:/nfs/software/helmod/apps/Core/R_core/3.6.1-gcb03/lib64/R/bin:/nfs/software/helmod/apps/Core/gcc/7.3.0-gcb01/bin:/nfs/software/helmod/apps/Core/pcre/8.37-fasrc02/bin:/nfs/software/helmod/apps/Core/pcre/8.37-fasrc02/bin:/nfs/software/helmod/apps/Core/xz/5.2.2-fasrc01/bin:/nfs/software/helmod/apps/Core/bzip2/1.0.6-fasrc01/bin:/nfs/software/helmod/apps/Core/curl/7.45.0-fasrc01/bin:/usr/lib64/qt-3.3/bin:/usr/local/bin:/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/usr/lpp/mmfs/bin:/opt/dell/srvadmin/bin:/home/asf40/bin:/usr/lpp/mmfs/bin:/opt/stack/bin:/data/tunglab/asf40/wgs_data/MedGenome_ftp/f4stats/AdmixTools-master/bin

# admixr tutorial: https://bodkan.net/admixr/articles/tutorial.html
# in R:
library(admixr, lib.loc="/data/tunglab/asf40/wgs_data/MedGenome_ftp/f4stats/AdmixTools-master/bin")
# set data prefix so R knows where to find our eigenstrat formatted files
data_prefix <- "./poop"
snps <- eigenstrat(data_prefix)
#snps <- eigenstrat(ind="/data/tunglab/asf40/wgs_data/MedGenome_ftp/f4stats/your.ind", snp="/data/tunglab/asf40/wgs_data/MedGenome_ftp/f4stats/all.chrX.snp", geno="/data/tunglab/asf40/wgs_data/MedGenome_ftp/f4stats/all.chrX.eigenstratgeno")
snps <- eigenstrat(ind="/data/tunglab/asf40/wgs_data/MedGenome_ftp/f4stats/scratch/tmp3.ind", snp="/data/tunglab/asf40/wgs_data/MedGenome_ftp/f4stats/scratch/tmp3.take2", geno="/data/tunglab/asf40/wgs_data/MedGenome_ftp/f4stats/scratch/tmp3.geno.original")

result <- d(W = c("Amboseli_historic", "Amboseli_recent"), X = "Mikumi", Y = "SW_olive", Z = "Macaque", data = snps)

read_ind(snps)
read_snp(snps, exclude = FALSE)
read_geno(snps)

logiinfo(result) # get run info


(prefix <- download_data(dirname = tempdir()))
list.files(path = dirname(prefix), pattern = basename(prefix), full.names = TRUE)
snps <- eigenstrat(prefix)
# in snps, we should have
#snps
#EIGENSTRAT object
#=================
#components:
#  ind file: /tmp/RtmpiHG45S/snps/snps.ind
#  snp file: /tmp/RtmpiHG45S/snps/snps.snp
#  geno file: /tmp/RtmpiHG45S/snps/snps.geno

# define populations
pops <- c("French", "Sardinian", "Han", "Papuan", "Khomani_San", "Mbuti", "Dinka")
# run f4ratio
result <- f4ratio(X = pops, A = "Altai", B = "Vindija", C = "Yoruba", O = "Chimp", data = snps)

# Rename HAP file
mv mapped.SRR2565914.wall.bam mapped.HAP.wall.bam

rm mapped.SRR256*.bam; rm mapped.SRR265*

# some fstats tutorials:
https://compvar-workshop.readthedocs.io/en/latest/contents/03_f3stats/f3stats.html

# plotting the results in R
tmp5 <- read.table("~/Downloads/f3_test", header=T)
tmp5$A_B_C <- paste(tmp5$A, tmp5$B, tmp5$C, sep="_")
# need to remove A-B/B-A duplicates
tmp6 <- tmp5[c(1:3)]
tmp7 <- tmp6[ !duplicated(apply(tmp6, 1, sort), MARGIN = 2), ]
tmp8 <- merge(tmp7, tmp5, by=c("A", "B", "C"), all.x=T)
ggplot()+geom_pointrange(data=tmp8, aes(x=reorder(A_B_C, f3), y=f3, ymin=f3-stderr, ymax=f3+stderr, color=Zscore), size=1)+theme_classic()+theme(text=element_text(size=14), axis.text.x = element_text(angle=90))+geom_hline(aes(yintercept=0), linetype="dashed", color="grey50")


```
