#! /bin/bash
#$ -wd /wynton/home/walllab/robinsonj
#$ -l h_rt=120:00:00,mem_free=7G,scratch=150G
#$ -pe smp 8
#$ -o /wynton/home/walllab/robinsonj/reports/baboon_preprocessing/WGSproc5
#$ -e /wynton/home/walllab/robinsonj/reports/baboon_preprocessing/WGSproc5


source /wynton/home/walllab/robinsonj/anaconda3/etc/profile.d/conda.sh
conda activate gentools2

# QSUB=/opt/sge/bin/lx-amd64/qsub


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

PROGRESSLOG=${SCRATCH}/${NAME}/${NAME}_WGSproc5_progress.log
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

cp ${SCRATCH}/${NAME}/${NAME}_BQSR.vcf.gz* .
exitVal=${?}

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    mv *log ${SCRATCH}/${NAME}
    exit
fi
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}


### Copy recalibration table
echo -e "[$(date "+%Y-%m-%d %T")] Copying recal.table1... " >> ${PROGRESSLOG}

cp ${SCRATCH}/${NAME}/${NAME}_BQSR_recal.table1 .
exitVal=${?}

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    mv *log ${SCRATCH}/${NAME}
    exit
fi
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}


### BQSR: make post-recalibration table
echo -e "[$(date "+%Y-%m-%d %T")] GATK BaseRecalibrator 2... " >> ${PROGRESSLOG}
LOG=${NAME}_H_d_recaltable2.log
date "+%Y-%m-%d %T" > ${LOG}

gatk3 -Xmx48g -Djava.io.tmpdir=./temp -T BaseRecalibrator \
-nct ${NSLOTS} \
-R ${REFERENCE} \
-I ${NAME}_realign.bam \
-knownSites ${NAME}_BQSR.vcf.gz \
-BQSR ${NAME}_BQSR_recal.table1 \
-o ${NAME}_BQSR_recal.table2 &>> ${LOG}

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    mv *log ${SCRATCH}/${NAME}
    exit
fi
date "+%Y-%m-%d %T" >> ${LOG}
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}


### BQSR: make recalibration plots
echo -e "[$(date "+%Y-%m-%d %T")] GATK AnalyzeCovariates... " >> ${PROGRESSLOG}
LOG=${NAME}_H_e_recalplots.log
date "+%Y-%m-%d %T" > ${LOG}

gatk3 -Xmx24g -Djava.io.tmpdir=./temp -T AnalyzeCovariates \
-R ${REFERENCE} \
-before ${NAME}_BQSR_recal.table1 \
-after ${NAME}_BQSR_recal.table2 \
-plots ${NAME}_BQSR_recal.plots.pdf

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    mv *log ${SCRATCH}/${NAME}
    mv *table2 ${SCRATCH}/${NAME}
    exit
fi
date "+%Y-%m-%d %T" >> ${LOG}
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}


### Save files to SCRATCH and remove realigned bam file
echo -e "[$(date "+%Y-%m-%d %T")] Saving recalibrated bam file... " >> ${PROGRESSLOG}

mv ${NAME}_BQSR_recal.plots.pdf ${SCRATCH}/${NAME}

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    mv *log ${SCRATCH}/${NAME}
    exit
fi
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}

mv *log ${SCRATCH}/${NAME}
rm ${SCRATCH}/${NAME}/${NAME}_BQSR.vcf.gz
rm ${SCRATCH}/${NAME}/${NAME}_BQSR.vcf.gz.tbi
rm ${SCRATCH}/${NAME}/${NAME}_BQSR_recal.table1
rm ${SCRATCH}/${NAME}/${NAME}_BQSR_recal.table2
rm ${SCRATCH}/${NAME}/${NAME}_realign.bam
rm ${SCRATCH}/${NAME}/${NAME}_realign.bai
rm ${SCRATCH}/${NAME}/${NAME}_realign.bam.intervals


### Clean up
cd ${TMPDIR}
rm -rf myworkdir
