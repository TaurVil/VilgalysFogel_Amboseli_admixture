#!/bin/bash

index=${SLURM_ARRAY_TASK_ID}
coverage=`head -$index 01_sets.txt | tail -1`
path_genome=./Reduced_Genome.fa

module load java; module load samtools; module load python
module load htslib
module load GATK
module load tabix

ls gVCF/*.$coverage*vcf.gz > 02.$coverage.list

GenomeAnalysisTK.sh CombineGVCFs  -R $path_genome -O test.$coverage.g.vcf.gz -V 02.$coverage.list
GenomeAnalysisTK.sh GenotypeGVCFs -R $path_genome -V test.$coverage.g.vcf.gz -O test.$coverage.vcf.gz --tmp-dir /data/tunglab/tpv/scratch/

rm test.$coverage.g.vcf.gz; rm test.$coverage.g.vcf.gz.tbi

##for chrom in `cut -f 1 01_targetted_chroms.bed`; do GenomeAnalysisTK.sh CombineGVCFs  -R $path_genome -O test.$coverage.$chrom.g.vcf.gz -V 02.$coverage.list -L $chrom; GenomeAnalysisTK.sh GenotypeGVCFs -R $path_genome -V test.$coverage.$chrom.g.vcf.gz -O test.$coverage.$chrom.vcf.gz --tmp-dir /data/tunglab/tpv/scratch/; rm test.$coverage.$chrom.g.vcf.gz; rm test.$coverage.$chrom.g.vcf.gz.tbi; done

tabix -f test.$coverage.vcf.gz
module load java/1.8.0_45-fasrc01; module load vcftools

java -jar ~/Programs/GenomeAnalysisTK.jar -T VariantFiltration -R $path_genome -V test.$coverage.vcf.gz -filterName "FS" --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" -o ./tmp.$coverage.vcf.gz

rm test.$coverage.vcf.gz
vcftools --gzvcf tmp.$coverage.vcf.gz --recode --out tmp.$coverage --max-alleles 2 --recode-INFO-all --remove-filtered-all

rm tmp.$coverage.vcf.gz; mv tmp.$coverage.recode.vcf test.$coverage.vcf; bgzip test.$coverage.vcf; tabix test.$coverage.vcf.gz

rm 02.$coverage.list;
