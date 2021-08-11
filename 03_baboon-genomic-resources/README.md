## Baboon genomic resources

This folder contains the baboon genomic resources (Panubis1 genome) which were used in this study, and the scripts used to generate them. 

* **Liftover files** were created between Panubis1 and hg38. Code to create the liftover files is available in `scripts/create-liftover-files`, which includes a separate README. Liftover files are available on Zenodo. 

* **Blood enhancers** were annotated by lifting over human PBMC H3K4me3 ChiP-seq peaks onto the baboon genome. Code used to do so is available in `scripts/Liftover_enhancers.sh` and the annotated peaks are in `resources/enhancers_lifted_to_baboon.bed`.

* **Recombination rates** were estimated across the genome using linkage information for 24 anubis SNPRC founder baboons using ldhelmet. The scripts to call this recombination map from raw genotype calls are included in `scripts/create-recombination-map`, which includes a separate README. This produces a file of recombination rates per chromosome, included in `resources/Recombination_Rates/`. Note that these are the population-scaled recombination rates returned by ldhelmet rather than a recombination rate in terms of morgans per base pair. 

* **B values** were estimated for the baboon genome using the recombination rate map called previously, a distribution of deleterious mutations estimated for humans (McVicker et al. 2009), and the location of genes in the baboon genome. Code to generate B values is available in `scripts/Annotate_B_values.sh` and the resulting files are available per chromsome in `resources`. There are two sets of output files, representing the decrease in neutral variation due to coding mutations (`B_values_coding`) and due to deleterious mutations in non-coding regions (`B_values_noncoding_near_genes`). 

* **Allele frequencies** for anubis and yellow baboons were estimated from all unadmixed individuals, after masking for putative heterospecific ancestry. This simply involved using `vcftools --freq` for sets of yellow and anubis baboons. Allele frequencies are available on Zenodo as `allele_frequencies.reference.masked.RData`. 

* **GC content** for 250 kb windows of the genome was calculated using bedtools. The resulting data file is included in `resources/gc_content.panubis1.250kb.txt`. 

* **CpG Islands** were annotated using EMBOSS with default settings. Code to do so from the Panubis1 genome is available in `scripts/Annotate_CpG_Islands.sh` and the resulting annotations are available in `resources/cpg_islands.panubis1.emboss.bed`. 