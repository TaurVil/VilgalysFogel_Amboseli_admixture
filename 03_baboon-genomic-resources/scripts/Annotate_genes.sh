## gene annotations performed by ncbi

#gtf-version 2.2
#!genome-build Panubis1.0
#!genome-build-accession NCBI_Assembly:GCF_008728515.1
#!annotation-source NCBI Papio anubis Annotation Release 104


## convert from NCBI RefSeq sequence names to the chromosomes: https://www.ncbi.nlm.nih.gov/assembly/GCF_008728515.1

sed -e 's/; / /g' Panubis1_chromnumbers.gtf > Panubis1_gtf_readable.txt
sed -i 's/gene_id "//g' Panubis1_gtf_readable.txt; sed -i 's/" transcript_id "/\t/g' Panubis1_gtf_readable.txt
sed -i 's/" db_xref "/\t/g' Panubis1_gtf_readable.txt; sed -i 's/" gbkey "/\t/g' Panubis1_gtf_readable.txt
sed -i 's/" gene "/\t/g' Panubis1_gtf_readable.txt; sed -i 's/" model_evidence "/\t/g' Panubis1_gtf_readable.txt
sed -i 's/" exon_number "/\t/g' Panubis1_gtf_readable.txt; sed -i 's/" product "/\t/g' Panubis1_gtf_readable.txt
sed -i 's/" protein_id "/\t/g' Panubis1_gtf_readable.txt; sed -i 's/" gene_biotype "/\t/g' Panubis1_gtf_readable.txt

sed -i 's/" exon_number "/\t/g' Panubis1_gtf_readable.txt; sed -i 's/" product "/\t/g' Panubis1_gtf_readable.txt