# Using code from McVicker et al. 2009

# Modify genotype calls/recombination map from original recombination output files 
	# So the file is recombination blocks, where the snp position is the start of this block and goes until the next block 
	# Need to make sure both the rate and length are in the right units 

	cp ./recombination/anubisSW.* ./

	module load R; R
    read.table("~/genomes/panubis1/Panubis1.0.fa.fai") -> lengths
	library(data.table) 
	# give lengths the length we'll be aiiming for in centiMorgans 
	# These are the chromosome lengths from Cox et al. 2006
	lengths$M <- NA; 
	lengths$M[1:20] <- c(172,115,100,164,113,111,116,93,88,88,91,126,77,76,69,89,77,84,63,63) 
	
	map <- NULL; for (i in 1:20) { 
		name <- paste("anubisSW.",i,".txt", sep="") 
        #name <- paste("output.",i,".txt", sep="") 
		d <- fread(name) 
		colnames(d)[1:6] <- c("left_snp", "right_snp", "mean", "lower", "median", "upper")
		d$Chromosome <- lengths$V1[i]
		d$`Position(bp)` <- d$left_snp #corrected to start, rather than midpoint:  round(rowMeans(d[,1:2])) # 
		
		# Remove recombination rates which are more than 100x the genome-wide average, setting them to that value 
		# Probably smoothes out the recombination map to some degree
        d$mean[d$mean > 100*7.71e-4] <- 100*7.71e-4
        print(c("Percent of windows smoothed", as.character(lengths$V1[i]), sum(d$mean == max(d$mean))/nrow(d)))
		## return percent of genome where RCR is limited 
		
		d$length_RCR <- (d$right_snp - d$left_snp)*d$mean 
		
		d2 <- subset(d, d$length_RCR > 0 )           #doesn't do anything now that the minimum length of a tract is 1bp  
		cumsum(d2$length_RCR) -> d2$map_pos1 # still in ldhelmet terms
		d2$map_pos1 * lengths$M[i] / max(d2$map_pos1) -> d2$`Map(cM)`
		d2$`Rate(cM/bp)` <- d2$mean * lengths$M[i] / max(d2$map_pos1)
		m2 <- rle(d2$`Rate(cM/bp)`) #rle pulls out the number of instances and number of times a value is repeated in a row. The cumsum of lengths (the number of repeats) is the unique instances. Works here since we're only defining where a tract starts
		Y <- cumsum(c(1, m2$lengths[-length(m2$lengths)]))
		d3 <- d2[Y,]

		name2 <- paste("panubis1_n24_genetic_map_",lengths$V1[i],".txt",sep="")
		write.table(d3[,c(8,9,13,12)], name2, row.names=F, col.names=T, sep="\t", quote=F)
		
		rbind(map, d2) -> map }
    
	## So I would need to adjust the output columns to use this... 
	# Scale to 2754M (total recombination map size estimated by Cox et al. 2006)
     map[,c(1:3,7:9,12,11),] -> map 
     cumsum(map$length_RCR) -> map$map
     map$`Rate2` <- map$mean * 2754 / max(map$map) #total recombination map length from Cox et al. 2006
    as.data.frame(map) -> map
    map2 <- NULL; for (i in 1:20) {
        d <- map[ map$Chromosome == as.character(lengths$V1[i]),] 
        d$length_RCR2 <- (d$right_snp - d$left_snp)*d$Rate2
        d$Map2 <- cumsum(d$length_RCR2)
        print(d[nrow(d),c(4,5,10,12,8,7)])
        m2 <- rle(d$`Rate2`)
        Y <- cumsum(c(1, m2$lengths[-length(m2$lengths)]))
        d3 <- d[Y,]
        rbind(map2, d3) -> map2 
        name2 <- paste("panubis1_n24_genetic_map2_",lengths$V1[i],".txt",sep="")
		write.table(d3[,c(4:5,10,12,8,7)], name2, row.names=F, col.names=T, sep="\t", quote=F); rm(m2,Y,d,d3) 
     }
    map2$contrib <- (map2$right_snp - map2$left_snp)*map2$mean
    map2$l <- (map2$right_snp - map2$left_snp)
    sum(map2$contrib)/sum(map2$l)
	write.table(map2, "../panubis1_n24_genetic_map_CoxTotalLength.txt", row.names=F, col.names=T, sep="\t", quote=F)

# get chromosome info file 
	cd /data/tunglab/tpv/B_values_GrahamCode/panubis1/; cp /data/tunglab/shared/genomes/panubis1/Panubis1.0.fa.fai ./chromInfo.txt 
	grep -v 'chrX' chromInfo.txt > tmp; grep -v 'chrY' tmp > chromInfo.txt; rm tmp 

# make separate conserved feature file for each chromosome 
	mkdir constrained_regions; cd constrained_regions; R
    library(data.table) 
	read.table("/data/tunglab/shared/genomes/panubis1/Panubis1.0.fa.fai") -> lengths
	fread("/data/tunglab/tpv/my_genomes/panubis1/Panubis1_chromnumbers.gtf", header=F) -> d
    d$V1 <- paste("chr", d$V1,sep=""); subset(d, d$V3 == 'exon') -> d
	for (i in 1:20) {
		subset(d, d$V1 == lengths$V1[i]) -> d2 
		d2[order(d2$V4,d2$V5),] -> d2 
		#out <- d %>%
		#	arrange(V2) %>%
		#	group_by(g = cumsum(cummax(lag(V3, default = first(V3))) < V2)) %>%
		#	summarise(start = first(V2), stop = max(V3))
        d2$V5[d2$V5 > lengths$V2[i]] <- lengths$V2[i]
        d2 <- subset(d2, d2$V5 - d2$V4 >=3)
		name=paste(lengths$V1[i], "_panubis1_exons.bed", sep="")
		write.table(d2[,-c(2:3,6:8)],name,col.names=F, row.names=F, sep="\t", quote=F)}
	fread("/data/tunglab/tpv/my_genomes/panubis1/Panubis1_chromnumbers.gtf", header=F) -> d
	d$V1 <- paste("chr", d$V1,sep=""); d <- subset(d, d$V3 == "gene")
    subset(d, d$V7 == '+') -> forw; subset(d, d$V7 == '-') -> rev
	rev$V5 <- rev$V5 + 1e4; forw$V4 <- forw$V4 - 1e4 #add 10kb promoter
	rbind(forw,rev) -> d; d[order(d$V1,d$V2),] -> d #combine forward and reverse again, sort 
	d$V4[d$V4 <= 0] <- 1
	for (i in 1:20) {
		subset(d, d$V1 == lengths$V1[i]) -> d2 
		d2[order(d2$V4,d2$V5),] -> d2 
		d2$V5[d2$V5 > lengths$V2[i]] <- lengths$V2[i]
		#out <- d %>%
		#	arrange(V2) %>%
		#	group_by(g = cumsum(cummax(lag(V3, default = first(V3))) < V2)) %>%
		#	summarise(start = first(V2), stop = max(V3))

		name=paste(lengths$V1[i], "_panubis1_genes_plus_10kbprom.bed", sep="")
		write.table(d2[,-c(2:3,6:8)],name,col.names=F, row.names=F, sep="\t", quote=F)
	}; quit(save='no')
	module load bedtools2; for f in `ls *bed`; do bedtools merge -i $f > tmp; mv tmp $f.v2; done
    cd /data/tunglab/tpv/B_values_GrahamCode/panubis1/


##
sbatch --mem=16G Annotate_B_values_script.sh


## Grab human data for comparison 
wget http://www.phrap.org/software_dir/mcvicker_dir/bkgd.tar.gz