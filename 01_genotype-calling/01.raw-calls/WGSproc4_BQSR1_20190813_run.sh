#! /bin/bash
#$ -wd /wynton/home/walllab/robinsonj
#$ -l h_rt=120:00:00,mem_free=7G,scratch=400G
#$ -pe smp 8
#$ -o /wynton/home/walllab/robinsonj/reports/baboon_preprocessing/WGSproc4
#$ -e /wynton/home/walllab/robinsonj/reports/baboon_preprocessing/WGSproc4


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

PROGRESSLOG=${SCRATCH}/${NAME}/${NAME}_WGSproc4_progress.log
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

cp ${SCRATCH}/${NAME}/${NAME}_realign.bam .
exitVal1=${?}
cp ${SCRATCH}/${NAME}/${NAME}_realign.bai .
exitVal2=${?}

if [ ${exitVal1} -ne 0 ] || [ ${exitVal2} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    exit
fi
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}


### Copy VCF files
echo -e "[$(date "+%Y-%m-%d %T")] Copying VCF files... " >> ${PROGRESSLOG}

cp ${SCRATCH}/${NAME}/${NAME}_BQSR_*.vcf.gz* .
exitVal=${?}

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    mv *log ${SCRATCH}/${NAME}
    exit
fi
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}


### Concatenate VCF files
echo -e "[$(date "+%Y-%m-%d %T")] Concatenating VCF files... " >> ${PROGRESSLOG}

ls ${NAME}_BQSR_*.vcf.gz > ${NAME}_BQSR_vcf.list
bcftools concat -f ${NAME}_BQSR_vcf.list -Oz -o ${NAME}_BQSR.vcf.gz
exitVal1=${?}
tabix -p vcf ${NAME}_BQSR.vcf.gz
exitVal2=${?}

if [ ${exitVal1} -ne 0 ] || [ ${exitVal2} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    exit
fi
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}

cp ${NAME}_BQSR.vcf.gz* ${SCRATCH}/${NAME}
rm ${NAME}_BQSR_*.vcf.gz*
rm ${SCRATCH}/${NAME}/${NAME}_BQSR_*.vcf.gz*


### BQSR: make recalibration table
echo -e "[$(date "+%Y-%m-%d %T")] GATK BaseRecalibrator 1... " >> ${PROGRESSLOG}
LOG=${NAME}_H_b_recaltable1.log
date "+%Y-%m-%d %T" > ${LOG}

gatk3 -Xmx48g -Djava.io.tmpdir=./temp -T BaseRecalibrator \
-nct ${NSLOTS} \
-R ${REFERENCE} \
-I ${NAME}_realign.bam \
-knownSites ${NAME}_BQSR.vcf.gz \
-o ${NAME}_BQSR_recal.table1 &>> ${LOG}

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    mv *log ${SCRATCH}/${NAME}
    exit
fi
date "+%Y-%m-%d %T" >> ${LOG}
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}

cp ${NAME}_BQSR_recal.table1 ${SCRATCH}/${NAME}


### BQSR: make recalibrated bam file
echo -e "[$(date "+%Y-%m-%d %T")] GATK PrintReads... " >> ${PROGRESSLOG}
LOG=${NAME}_H_c_printreads.log
date "+%Y-%m-%d %T" > ${LOG}

gatk3 -Xmx48g -Djava.io.tmpdir=./temp -T PrintReads \
-nct ${NSLOTS} \
-R ${REFERENCE} \
-BQSR ${NAME}_BQSR_recal.table1 \
--disable_indel_quals \
-I ${NAME}_realign.bam \
-o ${NAME}_recal.bam &>> ${LOG}

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    mv *log ${SCRATCH}/${NAME}
    exit
fi
date "+%Y-%m-%d %T" >> ${LOG}
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}


### Save files to SCRATCH
echo -e "[$(date "+%Y-%m-%d %T")] Saving recalibrated bam file... " >> ${PROGRESSLOG}

mv ${NAME}_recal.ba* ${SCRATCH}/${NAME}

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    mv *log ${SCRATCH}/${NAME}
    exit
fi
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}

mv *log ${SCRATCH}/${NAME}


### Submit next scripts
SCRIPTDIR=/wynton/home/walllab/robinsonj/scripts/WGSpipeline/baboon2019/single_highcov

SCRIPT=${SCRIPTDIR}/WGSproc5_BQSR2_20190813_run.sh
NEXT_JOB_ID1=$(${QSUB} -terse -N WGSproc5_${NAME} ${SCRIPT} ${NAME})
echo -e "[$(date "+%Y-%m-%d %T")] Submitted $(basename ${SCRIPT}): Job ID ${NEXT_JOB_ID1}" >> ${PROGRESSLOG}

SCRIPT=${SCRIPTDIR}/WGSproc6_HaplotypeCaller_20190813_run.sh
NJOBS=$(ls ${REFERENCE}.contiglist_* | wc -l)
NEXT_JOB_ID2=$(${QSUB} -terse -t 1-${NJOBS} -N WGSproc6_${NAME} ${SCRIPT} ${NAME} | cut -d'.' -f1)
echo -e "[$(date "+%Y-%m-%d %T")] Submitted $(basename ${SCRIPT}): Job ID ${NEXT_JOB_ID2}" >> ${PROGRESSLOG}

SCRIPT=${SCRIPTDIR}/WGSproc7_finalcheck_20190813_run.sh
NEXT_JOB_ID3=$(${QSUB} -terse -hold_jid ${NEXT_JOB_ID2} -N WGSproc7_${NAME} ${SCRIPT} ${NAME})
echo -e "[$(date "+%Y-%m-%d %T")] Submitted $(basename ${SCRIPT}) (holding on Job ID ${NEXT_JOB_ID2}): Job ID ${NEXT_JOB_ID3}" >> ${PROGRESSLOG}

SCRIPT=${SCRIPTDIR}/qualimap_20190813.sh
NEXT_JOB_ID4=$(${QSUB} -terse -N qmap_${NAME} ${SCRIPT} ${NAME})
echo -e "[$(date "+%Y-%m-%d %T")] Submitted $(basename ${SCRIPT}): Job ID ${NEXT_JOB_ID4}" >> ${PROGRESSLOG}


### Clean up
cd ${TMPDIR}
rm -rf myworkdir

