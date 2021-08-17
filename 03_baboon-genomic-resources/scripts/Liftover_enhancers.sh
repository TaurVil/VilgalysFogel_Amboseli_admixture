## To annotate enhancers in the baboon genome without direct experimental data, we'll liftover human enhancers to the baboon genome. We'll require 1:1 mapping such that each baboon enhancer maps back to the human genome in the same position which it originated from. Overall, we properly lift over 113k enhancers to the baboon genome, out of 188k human enhancers and 146k which initially lifted over to the baboon genome. 

######### Liftover enhancers
cd ~/genomes/humans
wget https://www.encodeproject.org/files/ENCFF995VPI/@@download/ENCFF995VPI.bed.gz
gunzip ENCFF995VPI.bed.gz; mv ENCFF995VPI.bed H3K4me1_peaks.bed

cd ~/genomes/panubis1
cut -f 1-4 ../humans/H3K4me1_peaks.bed > ./human_enhancer.bed

module load ucsc

## lift human H3K4me1_peaks over to the baboon genome (tmp.forward.bed)
liftOver ./human_enhancer.bed ~/my_genomes/liftover/hg38_to_Panubis1.0.chain.gz tmp.forward.bed forward_dropped_unMapped_enh.bed

## lift from baboon back to the human genome (tmp.reverse.bed)
liftOver ./tmp.forward.bed ~/my_genomes/liftover/Panubis1.0_to_hg38.chain.gz tmp.reverse.bed reverse_dropped_unMapped_enh.bed

## get only elements from tmp.forward.bed that are also in tmp.reverse.bed
awk 'NR==FNR{A[$1];next}$4 in A' tmp2.rev.names human_enhancer.bed > tmp3.bed 

## Restrict to elements which mapped to the same position 
module load R; R
read.delim("tmp3.bed", header=F) -> original_position; read.delim("tmp.forward.bed", header=F) -> forward_mapped; read.delim("tmp.reverse.bed", header=F) -> reverse_position
original_position$V4 <- as.character(original_position$V4); original_position$V1 <- as.character(original_position$V1)
exacts <- subset(original_position, original_position$V1 == reverse_position$V1 & original_position$V2 == reverse_position$V2 & original_position$V3 == reverse_position$V3 & original_position$V4 == reverse_position$V4)
f2 <- subset(forward_mapped, forward_mapped$V4 %in% exacts$V4)
write.table(f2, "./H3K4me1_peaks.bed", row.names=F, col.names=F, sep="\t", quote=F)

## on the github this file is renamed enhancers_lifted_to_baboon.bed 

## clean up
rm tmp3.bed; rm tmp3; rm tmp2.rev.names; rm tmp.forward.bed; rm tmp.reverse.bed; rm human_enhancer.bed; rm reverse_dropped_unMapped_enh.bed ; rm forward_dropped_unMapped_enh.bed

