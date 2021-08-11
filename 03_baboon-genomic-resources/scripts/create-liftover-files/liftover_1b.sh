#!/bin/bash
#$ -wd /data/tunglab/tpv/my_genomes/liftover/panubis1_to_hg38
#$ -l h_rt=24:00:00,mem_free=4G
#$ -N liftover1
#$ -V

# LiftOver pipeline, step 1: Preparing input fasta files and splitting

# Load relevant modules
module load ucsc; module load blat ; module load samtools


# USAGE:
# run this script

# Main directory location
MAINDIR=/data/tunglab/tpv/my_genomes/liftover/panubis1_to_hg38

# Specify old and new genomes
FASTADIR_NEW=/data/tunglab/shared/genomes/hg38_bwa/
NEW=hg38.fa

FASTADIR_OLD=/data/tunglab/shared/genomes/panubis1/
OLD=Panubis1.0.fa


# Generate .2bit and .sizes file for the OLD genome
mkdir -p ${MAINDIR}/fasta_old
cd ${MAINDIR}/fasta_old
cp ${FASTADIR_OLD}/${OLD} .
faToTwoBit ${OLD} ${OLD%.fa*}.2bit
twoBitInfo ${OLD%.fa*}.2bit ${OLD%.fa*}.sizes


# Generate .2bit and .sizes file for the NEW genome
mkdir -p ${MAINDIR}/fasta_new
cd ${MAINDIR}/fasta_new
cp ${FASTADIR_NEW}/${NEW} .
faToTwoBit ${NEW} ${NEW%.fa*}.2bit
twoBitInfo ${NEW%.fa*}.2bit ${NEW%.fa*}.sizes


# Split new genome by fasta record and generate .fai, .2bit, and .sizes for each part
mkdir -p ${MAINDIR}/fasta_new_split
cd ${MAINDIR}/fasta_new_split

faSplit byname ${FASTADIR_NEW}/${NEW} ./

for f in *.fa; do samtools faidx ${f} ; done
for f in *.fa; do faToTwoBit ${f%.*}.fa ${f%.*}.2bit ; done
for f in *.fa; do twoBitInfo ${f%.*}.2bit ${f%.*}.sizes ; done


# Generate .ooc file to speed up blat searches later
# Note: repMatch=1024 by default for tileSize=11, but should technically be set to
# repMatch = ( n / 2861349177 ) * 1024
# where n is equal to the number of "real" bases as reported by: faSize ${OLD}
cd ${MAINDIR}/fasta_old
tileSize=11
repMatch=1024
OOC=${OLD%.fa*}.${tileSize}.ooc
blat ${OLD%.fa*}.2bit /dev/null /dev/null -tileSize=${tileSize} -repMatch=${repMatch} -makeOoc=${OOC}


# Split individual fasta records into smaller chunks 
# (e.g. 200kb chunks for same species, 25kb chunks for different species)

# Same species:
# fastaSplitSize=200000

# Different species:
fastaSplitSize=25000

cd ${MAINDIR}
mkdir -p fasta_new_chunk

for f in ${MAINDIR}/fasta_new_split/*.fa ; do
    CHR=`basename ${f%.*}`
    DIR=${MAINDIR}/fasta_new_chunk/${CHR}_chunks
    mkdir -p ${DIR}
    cd ${DIR}
    faSplit size ${f} ${fastaSplitSize} ${CHR}_ -lift=${CHR}.lft
done

