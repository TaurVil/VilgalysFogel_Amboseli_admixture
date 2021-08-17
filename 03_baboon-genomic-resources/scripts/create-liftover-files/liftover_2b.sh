#!/bin/bash
#$ -l h_rt=24:00:00,mem_free=4G
#$ -N liftover2
#$ -V

# LiftOver pipeline, step 2: run blat on each chunk

# USAGE:

# Define variables from previous script
export MAINDIR=./panubis1_to_hg38
export FASTADIR_NEW=~/genomes/hg38_bwa/
export NEW=hg38.fa
export FASTADIR_OLD=~/genomes/panubis1
export OLD=Panubis1.0.fa
export tileSize=11
export OOC=${MAINDIR}/fasta_old/${OLD%.fa*}.${tileSize}.ooc

REPORTDIR=${MAINDIR}/reports; mkdir $REPORTDIR

# Submit blat job for each chunk   
for f in `ls ${MAINDIR}/fasta_new_split/*.fa` ; do
CHR=`basename ${f%.*}`
cd ${MAINDIR}/fasta_new_chunk/${CHR}_chunks
NUMJOBS=$(ls *.fa | wc -l)
sbatch --array=1-${NUMJOBS}%15 --mem=8G --nice ../../liftover_2b.subscript.sh
done

