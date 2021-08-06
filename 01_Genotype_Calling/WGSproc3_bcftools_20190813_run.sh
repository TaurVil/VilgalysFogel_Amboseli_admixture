#! /bin/bash
#$ -wd /wynton/home/walllab/robinsonj
#$ -l h_rt=48:00:00,mem_free=24G,scratch=150G
#$ -o /wynton/home/walllab/robinsonj/reports/baboon_preprocessing/WGSproc3
#$ -e /wynton/home/walllab/robinsonj/reports/baboon_preprocessing/WGSproc3


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

PROGRESSLOG=${SCRATCH}/${NAME}/${NAME}_WGSproc3_${IDX}_progress.log
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


### Generate raw VCF file for use in recalibration
echo -e "[$(date "+%Y-%m-%d %T")] bcftools mpileup/call/filter... " >> ${PROGRESSLOG}
LOG=${NAME}_H_a_bcftools_${IDX}.log
date "+%Y-%m-%d %T" > ${LOG}

set -o pipefail

bcftools mpileup -f ${REFERENCE} -R ${REFERENCE}.contiglist_${IDX}.list --annotate INFO/ADF,INFO/ADR \
-Q 20 -q 30 --rf 2 --ff 256 --ff 1024 -Ou ${NAME}_realign.bam 2>> ${LOG} | \
bcftools call -mv -Ou 2>> ${LOG} | \
bcftools filter -i 'N_ALT=1 && INFO/ADF[1]>0 && INFO/ADR[1]>0' \
-Oz -o ${NAME}_BQSR_${IDX}.vcf.gz 2>> ${LOG}

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    mv *log ${SCRATCH}/${NAME}
    exit
fi
date "+%Y-%m-%d %T" >> ${LOG}
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}


### Index VCF file
echo -e "[$(date "+%Y-%m-%d %T")] tabix indexing... " >> ${PROGRESSLOG}

tabix -p vcf ${NAME}_BQSR_${IDX}.vcf.gz

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    mv *log ${SCRATCH}/${NAME}
    exit
fi
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}


### Save files to SCRATCH
echo -e "[$(date "+%Y-%m-%d %T")] Saving VCF file for BQSR... " >> ${PROGRESSLOG}

mv ${NAME}_BQSR_${IDX}.vcf.gz* ${SCRATCH}/${NAME}

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
