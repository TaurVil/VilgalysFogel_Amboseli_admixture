**ADLIBS (Schaefer et al. 2017)** uses a Hidden Markov Model (HMM) and pseudo-haploid chromosomes to identify regions of the genome that more closely resemble one parental population than the other. Because this method uses pseudo-haploid data, it can be inferred from a single base pair at a position and additional coverage doesn't affect the outcome. However, ADLIBS appears to perform poorly when faced with abundant segregating variation, because haploid information is unable to reliably assign ancestry in the absence of fixed or nearly fixed differences. 

```console
##haploidize each indivdiual
###Helpful hints for running this from: https://github.com/jacahill/Admixture
###For this pipeline you should have one fully filtered bam file per sample. haploidize each bam file to a fasta. 
###I use R. Ed Green's pu2fa program available here- https://github.com/Paleogenomics/Chrom-Compare Usage: For each chromosome run:

## for test individuals at different coverage
for f in `cat 00names`; do for g in `cat 01_sets.txt`; do for h in `cat 00chroms`; do cat get_haploid_fasta.sh | sed -e s/SCAF/$h/g | sed -e s/NAME/$f/g | sed -e s/COV/$g/g > g.sh; sbatch --mem=13000 g.sh; done; done; done

## for reference panel individuals
for f in `cat 00refpanel`; do for h in `cat 00chroms`; do cat get_haploid_fasta_reference_.sh | sed -e s/SCAF/$h/g | sed -e s/NAME/$f/g  > g.sh; sbatch --mem=13000 g.sh; done; done

## Run ADLIBS and then extract ancestry tracts
for f in `cat 01_sets.txt`; do for g in `cat 00chroms`; do cat run_adlibs.sh | sed -e s/CHROM/$g/g | sed -e s/XXX/$f/g > g.sh; sbatch --nice --mem=13000 g.sh; done; done 

## for each coverage and individual, merge together ancestry tracts per chromosome
for f in `cat 00covs`; do for g in `cat 00names`; do cat run1.$f.*/$f.$g_*.bed > $g.$f.25k.50.bed; done; done

```
