#!/bin/bash

NAME=${1}

FASTQDIR=/mnt/Wall1/baboon/fastq/${NAME}

WYNTONSCRIPT=/wynton/home/walllab/robinsonj/scripts/WGSpipeline/baboon2019/single_highcov/WGSproc1_trim_align_sort_addRG_20190813_run.sh

FQ1=$(ls ${FASTQDIR}/*gz | head -n 1)
FQ2=$(ls ${FASTQDIR}/*gz | tail -n 1)
R1=$(basename ${FQ1})
R2=$(basename ${FQ2})
RGID=${NAME}_01
RGLB=Lib01
RGPU=$(zcat ${FQ1} | head -n 1 | cut -d':' -f3-4 | sed 's/:/./g' ).${NAME}
RGPL=illumina
FLAG=1

ssh robinsonj@log2.wynton.ucsf.edu "rm -rf /wynton/scratch/robinsonj/baboonWGS/${NAME}"
ssh robinsonj@log2.wynton.ucsf.edu "mkdir -p /wynton/scratch/robinsonj/baboonWGS/${NAME}"

scp ${FQ1} robinsonj@dt1.wynton.ucsf.edu:/wynton/scratch/robinsonj/baboonWGS/${NAME}
scp ${FQ2} robinsonj@dt1.wynton.ucsf.edu:/wynton/scratch/robinsonj/baboonWGS/${NAME}

echo -e "[$(date "+%Y-%m-%d %T")] qsub ${WYNTONSCRIPT} ${R1} ${R2} ${NAME} ${RGID} ${RGLB} ${RGPU} ${RGPL} ${FLAG}"

ssh robinsonj@log2.wynton.ucsf.edu "qsub -N WGSproc1_${NAME} ${WYNTONSCRIPT} ${R1} ${R2} ${NAME} ${RGID} ${RGLB} ${RGPU} ${RGPL} ${FLAG}"
