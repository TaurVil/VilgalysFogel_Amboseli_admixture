#!/bin/bash
#$ -wd /data/tunglab/tpv/my_genomes/liftover/hg38_to_panubis1
#$ -l h_rt=24:00:00,mem_free=4G
#$ -N liftover2
#$ -V

# LiftOver pipeline, step 2: run blat on each chunk

# USAGE:
# run this script 


# Define variables from previous script
export MAINDIR=/data/tunglab/tpv/scratch/hg38_to_panubis1_liftover
export FASTADIR_OLD=/data/tunglab/shared/genomes/hg38_bwa/
export OLD=hg38.fa
export FASTADIR_NEW=/data/tunglab/shared/genomes/panubis1
export NEW=Panubis1.0.fa
export tileSize=11
export OOC=${MAINDIR}/fasta_old/${OLD%.fa*}.${tileSize}.ooc

REPORTDIR=${MAINDIR}/reports; mkdir $REPORTDIR

# Submit blat job for each chunk   ---- this is the original code from Jacqueline which we need to modify for hardac
for f in `ls ${MAINDIR}/fasta_new_split/*.fa` ; do
CHR=`basename ${f%.*}`
cd ${MAINDIR}/fasta_new_chunk/${CHR}_chunks
NUMJOBS=$(ls *.fa | wc -l)
sbatch --array=1-${NUMJOBS}%45 --mem=8G ../../liftover_2.subscript.sh
done



