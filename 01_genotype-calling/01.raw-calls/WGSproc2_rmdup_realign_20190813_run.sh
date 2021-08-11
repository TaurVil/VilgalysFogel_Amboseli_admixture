#! /bin/bash
#$ -wd /wynton/home/walllab/robinsonj
#$ -l h_rt=96:00:00,mem_free=56G,scratch=350G
#$ -o /wynton/home/walllab/robinsonj/reports/baboon_preprocessing/WGSproc2
#$ -e /wynton/home/walllab/robinsonj/reports/baboon_preprocessing/WGSproc2


source /wynton/home/walllab/robinsonj/anaconda3/etc/profile.d/conda.sh
conda activate gentools2

QSUB=/opt/sge/bin/lx-amd64/qsub


### Variables
NAME=${1}

SCRATCH=/wynton/scratch/robinsonj/baboonWGS

REFDIR=/wynton/home/walllab/robinsonj/work/baboon/reference/Panubis_1.0
REFERENCE=Panubis_1.0.fasta


### Use TMPDIR as the working directory
cd ${TMPDIR}
mkdir myworkdir
cd myworkdir
mkdir temp

PROGRESSLOG=${SCRATCH}/${NAME}/${NAME}_WGSproc2_progress.log
echo -e "[$(date "+%Y-%m-%d %T")] JOB_ID: ${JOB_ID}" > ${PROGRESSLOG}


### Copy reference genome
echo -e "[$(date "+%Y-%m-%d %T")] Copying reference genome... " >> ${PROGRESSLOG}

cp ${REFDIR}/* .

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    exit
fi
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}


### Copy alignment files
echo -e "[$(date "+%Y-%m-%d %T")] Copying alignment files... " >> ${PROGRESSLOG}

cp ${SCRATCH}/${NAME}/${NAME}_*_addRG.bam .
exitVal1=${?}
cp ${SCRATCH}/${NAME}/${NAME}_*_addRG.bai .
exitVal2=${?}

if [ ${exitVal1} -ne 0 ] || [ ${exitVal2} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    exit
fi
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}


### Remove duplicates
echo -e "[$(date "+%Y-%m-%d %T")] Picard MarkDuplicates... " >> ${PROGRESSLOG}
LOG=${NAME}_F_rmDup.log
date "+%Y-%m-%d %T" > ${LOG}

picard -Xmx48G MarkDuplicates \
$(for i in *addRG.bam; do echo "INPUT=${i} "; done) \
OUTPUT=${NAME}_rmDup.bam \
METRICS_FILE=${NAME}_rmDup_metrics.txt \
MAX_RECORDS_IN_RAM=250000 MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
REMOVE_DUPLICATES=true \
CREATE_INDEX=true TMP_DIR=./temp 2>> ${LOG}

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    mv *log ${SCRATCH}/${NAME}
    exit
fi
date "+%Y-%m-%d %T" >> ${LOG}
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}

mv ${NAME}_rmDup_metrics.txt ${SCRATCH}/${NAME}
rm ${NAME}_*_addRG.ba*


### Indel realignment A: create target intervals
echo -e "[$(date "+%Y-%m-%d %T")] GATK RealignerTargetCreator... " >> ${PROGRESSLOG}
LOG=${NAME}_G_a_realign.log
date "+%Y-%m-%d %T" > ${LOG}

gatk3 -Xmx48G -Djava.io.tmpdir=./temp -T RealignerTargetCreator \
-nt ${NSLOTS} \
-R ${REFERENCE} \
-I ${NAME}_rmDup.bam \
-o ${NAME}_realign.bam.intervals 2>> ${LOG}

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    mv *log ${SCRATCH}/${NAME}
    exit
fi
date "+%Y-%m-%d %T" >> ${LOG}
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}


### Indel realignment B: realign
echo -e "[$(date "+%Y-%m-%d %T")] GATK IndelRealigner... " >> ${PROGRESSLOG}
LOG=${NAME}_G_b_realign.log
date "+%Y-%m-%d %T" > ${LOG}

gatk3 -Djava.io.tmpdir=./temp -Xmx48G -T IndelRealigner \
-R ${REFERENCE} \
-I ${NAME}_rmDup.bam \
-o ${NAME}_realign.bam \
-targetIntervals ${NAME}_realign.bam.intervals 2>> ${LOG}

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    mv *log ${SCRATCH}/${NAME}
    exit
fi
date "+%Y-%m-%d %T" >> ${LOG}
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}

rm ${NAME}_rmDup.ba*


### Save files to SCRATCH and remove raw alignments
echo -e "[$(date "+%Y-%m-%d %T")] Saving realigned bam file... " >> ${PROGRESSLOG}

mv ${NAME}_realign.ba* ${SCRATCH}/${NAME}

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    mv *log ${SCRATCH}/${NAME}
    exit
fi
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}

mv *log ${SCRATCH}/${NAME}
rm ${SCRATCH}/${NAME}/${NAME}_*_addRG.bam
rm ${SCRATCH}/${NAME}/${NAME}_*_addRG.bai


### Submit next scripts
SCRIPTDIR=/wynton/home/walllab/robinsonj/scripts/WGSpipeline/baboon2019/single_highcov

SCRIPT=${SCRIPTDIR}/WGSproc3_bcftools_20190813_run.sh
NJOBS=$(ls ${REFERENCE}.contiglist_* | wc -l)
NEXT_JOB_ID1=$(${QSUB} -terse -t 1-${NJOBS} -N WGSproc3_${NAME} ${SCRIPT} ${NAME} | cut -d'.' -f1)
echo -e "[$(date "+%Y-%m-%d %T")] Submitted $(basename ${SCRIPT}): Job ID ${NEXT_JOB_ID1}" >> ${PROGRESSLOG}

SCRIPT=${SCRIPTDIR}/WGSproc4_BQSR1_20190813_run.sh
NEXT_JOB_ID2=$(${QSUB} -terse -hold_jid ${NEXT_JOB_ID1} -N WGSproc4_${NAME} ${SCRIPT} ${NAME})
echo -e "[$(date "+%Y-%m-%d %T")] Submitted $(basename ${SCRIPT}) (holding on Job ID ${NEXT_JOB_ID1}): Job ID ${NEXT_JOB_ID2}" >> ${PROGRESSLOG}


### Clean up
cd ${TMPDIR}
rm -rf myworkdir
