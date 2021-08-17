#!/bin/bash
#$ -l h_rt=24:00:00,mem_free=4G
#$ -N liftover4
#$ -V

# LiftOver pipeline, step 4: merging chain files, netting, and creating final chain

# USAGE:
# SCRIPT=./liftover_4.sh
# sbatch ${SCRIPT}

#HARDAC locations
MAINDIR=./panubis1_to_hg38
FASTADIR_NEW=~/genomes/hg38_bwa/
NEW=hg38.fa
FASTADIR_OLD=~/genomes/panubis1
OLD=Panubis1.0.fa

#Load HARDAC MODULES
module load ucsc; module load samtools

# Merge individual sorted chains
cd ${MAINDIR}/chains
chainMergeSort per_fasta_sorted/*sorted.chain | chainSplit . stdin -lump=1
mv 000.chain ${OLD%.fa}_to_${NEW%.fa}_merged.chain


# Make net file
chainNet ${OLD%.fa}_to_${NEW%.fa}_merged.chain ${MAINDIR}/fasta_old/${OLD%.fa*}.sizes ${MAINDIR}/fasta_new/${NEW%.fa*}.sizes ${OLD%.fa}_to_${NEW%.fa}.net /dev/null


# Make final chain file
netChainSubset ${OLD%.fa}_to_${NEW%.fa}.net ${OLD%.fa}_to_${NEW%.fa}_merged.chain ${OLD%.fa}_to_${NEW%.fa}.chain


# Compress final chain file (gzip)
gzip ${OLD%.fa}_to_${NEW%.fa}.chain
mv ${OLD%.fa}_to_${NEW%.fa}.chain.gz ${MAINDIR}/${OLD%.fa}_to_${NEW%.fa}.chain.gz


