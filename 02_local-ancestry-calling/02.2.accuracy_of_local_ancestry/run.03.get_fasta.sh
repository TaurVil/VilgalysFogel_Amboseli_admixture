#!/bin/bash

module load python/2.7.6-fasrc01;

## add header to the vcf file
cat full_header.vcf simulated_vcfs/INDIV.CHROMO.vcf | grep -v '0\/0' > simulated_vcfs/fixed.INDIV.CHROMO.vcf
sed -i 's/CHROM/#CHROM/g' simulated_vcfs/fixed.INDIV.CHROMO.vcf ; sed -i 's/\tQD/\tPASS/g' simulated_vcfs/fixed.INDIV.CHROMO.vcf

python ~/Programs/neat-genreads/genReads.py -r ./partial_genome/CHROMO.fa -R 101 -o ./sim_reads/INDIV.CHROMO -c 10 -v ./simulated_vcfs/fixed.INDIV.CHROMO.vcf -M 0 --rng 10 --gz --vcf

#could also merge vcfs and use Reduced_Genome.fa, but it will be faster to do by chromosome

# alternate options and parameter explanations
#c adjusts the average coverage from a 10x default
#leaving ploidy at 2
#using VCF to insert mutations rather than generating them. -v is the mutations to insert. -M is the mutation rate beyond those.
#--rng is the set-seed value
#paired end (pe) data with mean length 300 and sd of 30 could be used instead of se (`-R LENGTH`) by using: --pe 300 30
#output true vcf: --vcf

