#! /bin/bash
#$ -wd /wynton/home/walllab/robinsonj
#$ -l h_rt=120:00:00,mem_free=8G,scratch=500G
#$ -pe smp 4
#$ -o /wynton/home/walllab/robinsonj/reports/baboon_preprocessing/WGSproc1
#$ -e /wynton/home/walllab/robinsonj/reports/baboon_preprocessing/WGSproc1


source /wynton/home/walllab/robinsonj/anaconda3/etc/profile.d/conda.sh
conda activate gentools2

QSUB=/opt/sge/bin/lx-amd64/qsub


### Variables
FQ1=${1}
FQ2=${2}
NAME=${3}
RGID=${4}
RGLB=${5}
RGPU=${6}
RGPL=${7}
FLAG=${8}

SCRATCH=/wynton/scratch/robinsonj/baboonWGS

REFDIR=/wynton/home/walllab/robinsonj/work/baboon/reference/Panubis_1.0
REFERENCE=Panubis_1.0.fasta


### Use TMPDIR as the working directory
cd ${TMPDIR}
mkdir myworkdir
cd myworkdir
mkdir temp

PROGRESSLOG=${SCRATCH}/${NAME}/${RGID}_WGSproc1_progress.log
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


### Copy reads
echo -e "[$(date "+%Y-%m-%d %T")] Copying fastq files... " >> ${PROGRESSLOG}

cp ${SCRATCH}/${NAME}/${FQ1} .
exitVal1=${?}
cp ${SCRATCH}/${NAME}/${FQ2} .
exitVal2=${?}

if [ ${exitVal1} -ne 0 ] || [ ${exitVal2} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    exit
fi
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}


### FastQC for raw reads
echo -e "[$(date "+%Y-%m-%d %T")] Running FastQC on raw reads... " >> ${PROGRESSLOG}

fastqc -t ${NSLOTS} ${FQ1} ${FQ2}

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    exit
fi
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}


### Trim adapter sequences
echo -e "[$(date "+%Y-%m-%d %T")] BBDUK... " >> ${PROGRESSLOG}
LOG=${RGID}_A_bbduk.log
date "+%Y-%m-%d %T" > ${LOG}

bbduk.sh -Xmx26G threads=${NSLOTS} ref=${BBDUK_REF} \
ktrim=r k=23 mink=11 hdist=1 tpe tbo mlf=0.5 \
in1=${FQ1} in2=${FQ2} \
out1=${RGID}_R1_trim.fastq.gz out2=${RGID}_R2_trim.fastq.gz \
stats=${RGID}_A_bbduk_report.txt 2>> ${LOG}

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    mv *log ${SCRATCH}/${NAME}
    exit
fi
date "+%Y-%m-%d %T" >> ${LOG}
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}


### FastQC for trimmed reads
echo -e "[$(date "+%Y-%m-%d %T")] Running FastQC on trimmed reads... " >> ${PROGRESSLOG}

fastqc -t ${NSLOTS} ${RGID}_R1_trim.fastq.gz ${RGID}_R2_trim.fastq.gz

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    exit
fi
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}

mv *fastqc* ${SCRATCH}/${NAME}
mv ${RGID}_A_bbduk_report.txt ${SCRATCH}/${NAME}
rm ${FQ1}
rm ${FQ2}


### Align trimmed reads to reference and sort by query name
echo -e "[$(date "+%Y-%m-%d %T")] BWA MEM/samtools sort by query... " >> ${PROGRESSLOG}
LOG1=${RGID}_B_bwamem.log
date "+%Y-%m-%d %T" > ${LOG1}
LOG2=${RGID}_C_querysort.log
date "+%Y-%m-%d %T" > ${LOG2}

READ1=${RGID}_R1_trim.fastq.gz
READ2=${RGID}_R2_trim.fastq.gz

set -o pipefail

bwa mem -M -t ${NSLOTS} ${REFERENCE} ${READ1} ${READ2} 2>> ${LOG1} | \
samtools sort -n -@ ${NSLOTS} -O BAM -T ./temp/${NAME}_${RGID} \
-o ${RGID}_querysort.bam - 2>> ${LOG2}

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    mv *log ${SCRATCH}/${NAME}
    exit
fi
date "+%Y-%m-%d %T" >> ${LOG1}
date "+%Y-%m-%d %T" >> ${LOG2}
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}

rm ${READ1}
rm ${READ2}


### Fix mate information
echo -e "[$(date "+%Y-%m-%d %T")] samtools fixmate... " >> ${PROGRESSLOG}
LOG=${RGID}_D_fixmate.log
date "+%Y-%m-%d %T" > ${LOG}

samtools fixmate -@ ${NSLOTS} -O BAM ${RGID}_querysort.bam \
${RGID}_fixmate.bam 2>> ${LOG}

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    mv *log ${SCRATCH}/${NAME}
    exit
fi
date "+%Y-%m-%d %T" >> ${LOG}
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}

rm ${RGID}_querysort.bam


### Sort by coordinate and add read group info
echo -e "[$(date "+%Y-%m-%d %T")] Picard AddOrReplaceReadGroups... " >> ${PROGRESSLOG}
LOG=${RGID}_E_addRG.log
date "+%Y-%m-%d %T" > ${LOG}

picard -Xmx26G AddOrReplaceReadGroups \
INPUT=${RGID}_fixmate.bam OUTPUT=${RGID}_addRG.bam \
RGID=${RGID} RGSM=${NAME} RGLB=${RGLB} RGPU=${RGPU} RGPL=${RGPL} \
SORT_ORDER=coordinate CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT \
TMP_DIR=./temp 2>> ${LOG}

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    mv *log ${SCRATCH}/${NAME}
    exit
fi
date "+%Y-%m-%d %T" >> ${LOG}
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}

rm ${RGID}_fixmate.bam


### Save files to SCRATCH and remove raw reads
echo -e "[$(date "+%Y-%m-%d %T")] Saving addRG bam file... " >> ${PROGRESSLOG}

mv ${RGID}_addRG.ba* ${SCRATCH}/${NAME}

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    mv *log ${SCRATCH}/${NAME}
    exit
fi
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}

mv *log ${SCRATCH}/${NAME}
rm ${SCRATCH}/${NAME}/${FQ1}
rm ${SCRATCH}/${NAME}/${FQ2}


### Submit next script if only one read set for sample (flag=0: do not qsub, flag=1: qsub)
if [ ${FLAG} -ne 0 ]; then
    SCRIPTDIR=/wynton/home/walllab/robinsonj/scripts/WGSpipeline/baboon2019/single_highcov
    SCRIPT=${SCRIPTDIR}/WGSproc2_rmdup_realign_20190813_run.sh
    NEXT_JOB_ID=$(${QSUB} -terse -N WGSproc2_${NAME} ${SCRIPT} ${NAME})
    echo -e "[$(date "+%Y-%m-%d %T")] Submitted $(basename ${SCRIPT}): Job ID ${NEXT_JOB_ID}" >> ${PROGRESSLOG}
fi


### Clean up
cd ${TMPDIR}
rm -rf myworkdir
