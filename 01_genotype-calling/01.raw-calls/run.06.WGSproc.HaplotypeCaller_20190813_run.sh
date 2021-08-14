#! /bin/bash
#$ -wd /wynton/home/walllab/robinsonj
#$ -l h_rt=72:00:00,mem_free=48G,scratch=150G
#$ -o /wynton/home/walllab/robinsonj/reports/baboon_preprocessing/WGSproc6
#$ -e /wynton/home/walllab/robinsonj/reports/baboon_preprocessing/WGSproc6


source /wynton/home/walllab/robinsonj/anaconda3/etc/profile.d/conda.sh
conda activate gentools2

# QSUB=/opt/sge/bin/lx-amd64/qsub


### Variables
NAME=${1}
IDX=$(printf %02d ${SGE_TASK_ID})

SCRATCH=/wynton/scratch/robinsonj/baboonWGS

REFDIR=/wynton/home/walllab/robinsonj/work/baboon/reference/Panubis_1.0
REFERENCE=Panubis_1.0.fasta


### Use TMPDIR as the working directory
cd ${TMPDIR}
mkdir myworkdir
cd myworkdir
mkdir temp

PROGRESSLOG=${SCRATCH}/${NAME}/${NAME}_WGSproc6_${IDX}_progress.log
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

cp ${SCRATCH}/${NAME}/${NAME}_recal.bam .
exitVal1=${?}
cp ${SCRATCH}/${NAME}/${NAME}_recal.bai .
exitVal2=${?}

if [ ${exitVal1} -ne 0 ] || [ ${exitVal2} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    exit
fi
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}


### Generate gVCF files
echo -e "[$(date "+%Y-%m-%d %T")] GATK HaplotypeCaller... " >> ${PROGRESSLOG}
LOG=${NAME}_I_HaplotypeCaller_${IDX}.log
date "+%Y-%m-%d %T" > ${LOG}

gatk3 -Xmx40g -Djava.io.tmpdir=./temp -T HaplotypeCaller \
-nct ${NSLOTS} \
-R ${REFERENCE} \
-ERC GVCF \
-mbq 20 \
-I ${NAME}_recal.bam \
-L ${REFERENCE}.contiglist_${IDX}.list \
-o ${NAME}_${IDX}.g.vcf.gz &>> ${LOG}

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    mv *log ${SCRATCH}/${NAME}
    exit
fi
date "+%Y-%m-%d %T" >> ${LOG}
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}


### Save files to SCRATCH
echo -e "[$(date "+%Y-%m-%d %T")] Saving gVCF file... " >> ${PROGRESSLOG}

mv ${NAME}_${IDX}.g.vcf.gz* ${SCRATCH}/${NAME}

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    mv *log ${SCRATCH}/${NAME}
    exit
fi
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}

mv *log ${SCRATCH}/${NAME}


### Clean up
cd ${TMPDIR}
rm -rf myworkdir

