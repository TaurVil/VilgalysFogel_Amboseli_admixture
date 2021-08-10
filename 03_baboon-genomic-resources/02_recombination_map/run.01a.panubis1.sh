#!/bin/bash

module load vcftools; module load tabix; module load samtools
samtools faidx ../Panubis1.nochromname.fa CHROM | vcf-consensus -H 1 -s INDIV CHROM.anubisSW.recode.vcf.gz > individual_haplotype_fastas/INDIV.CHROM.1.fa
samtools faidx ../Panubis1.nochromname.fa CHROM | vcf-consensus -H 2 -s INDIV CHROM.anubisSW.recode.vcf.gz > individual_haplotype_fastas/INDIV.CHROM.2.fa

sed -i 's/>CHROM/>CHROM_INDIV.1/g' individual_haplotype_fastas/INDIV.CHROM.1.fa
sed -i 's/>CHROM/>CHROM_INDIV.2/g' individual_haplotype_fastas/INDIV.CHROM.2.fa
