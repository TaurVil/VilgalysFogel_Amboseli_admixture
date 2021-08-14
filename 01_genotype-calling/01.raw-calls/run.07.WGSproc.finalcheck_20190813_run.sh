#! /bin/bash
#$ -wd /wynton/home/walllab/robinsonj
#$ -l h_rt=00:05:00,mem_free=1G
#$ -o /wynton/home/walllab/robinsonj/reports/baboon_preprocessing/WGSproc7
#$ -e /wynton/home/walllab/robinsonj/reports/baboon_preprocessing/WGSproc7


### Variables
NAME=${1}

SCRATCH=/wynton/scratch/robinsonj/baboonWGS
SAVEDIR=/wynton/home/walllab/robinsonj/work/baboon/WGSpreprocessing


### Check to make sure that there were no reported failures
cd ${SCRATCH}/${NAME}
FCOUNT=$(cat *progress.log | grep "FAIL" | wc -l)

if [ $FCOUNT -eq 0 ]
then
    echo ${NAME} >> ${SAVEDIR}/donelist_single_highcov.txt
else
    echo ${NAME} >> ${SAVEDIR}/faillist_single_highcov.txt
fi

