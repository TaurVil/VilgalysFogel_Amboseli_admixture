#!/bin/bash
#$ -l h_rt=24:00:00,mem_free=4G
#$ -N liftover3
#$ -V

# LiftOver pipeline, step 3: creating final chain files from fragmented .psl files

# USAGE -- HARDAC
# export MAINDIR=/data/tunglab/tpv/my_genomes/liftover/hg38_to_panubis1
# cd ${MAINDIR}/fasta_new_split
# NUMJOBS=$(ls *.fa | wc -l)
# SCRIPT=/data/tunglab/tpv/my_genomes/liftover/hg38_to_panubis1/liftover_3.sh
# sbatch -a 1-${NUMJOBS} --mem=28G ${SCRIPT}

#Pulling these from previous scripts --HARDAC
MAINDIR=/data/tunglab/tpv/my_genomes/liftover/hg38_to_panubis1
FASTADIR_OLD=/data/tunglab/shared/genomes/hg38_bwa/
OLD=hg38.fa
FASTADIR_NEW=/data/tunglab/shared/genomes/panubis1
NEW=Panubis1.0.fa

#Get chromosome/scaffold name-- HARDAC
cd ${MAINDIR}/fasta_new_split
CHR=$(ls *.fa | head -${SLURM_ARRAY_TASK_ID} | tail -1 | sed 's/.fa//')

#Load necessary modules
module load ucsc

# Lift and merge per-chunk .psl files
cd ${MAINDIR}/fasta_new_chunk/${CHR}_chunks
liftUp -type=.psl -pslQ ${CHR}_lifted.psl ${CHR}.lft warn *temp.psl


# Make chain files
axtChain -linearGap=medium -psl ${CHR}_lifted.psl ${MAINDIR}/fasta_old/${OLD%.fa*}.2bit ${MAINDIR}/fasta_new_split/${CHR}.2bit ${CHR}.chain
chainSort ${CHR}.chain ${CHR}_sorted.chain


# Make chain output folder if it does not already exist and move sorted chain file
mkdir -p ${MAINDIR}/chains
mkdir -p ${MAINDIR}/chains/per_fasta_sorted
mv ${CHR}_sorted.chain ${MAINDIR}/chains/per_fasta_sorted

