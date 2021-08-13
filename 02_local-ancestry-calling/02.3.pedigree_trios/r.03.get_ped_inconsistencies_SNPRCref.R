#!/bin/env Rscript

# Load R data file for the masked SNPRC reference panel generated in 1prep_for_pedigree_inconsistencies.R script
load("local_ancestry_pedigree_trios_maskedSNPRCref.Rd")

# Load all ancestry calls for all individuals in the pedigree trios at all sites of interest across the genome
calls <- read.table("all.ancestry_calls_maskedSNPRCref_pedigree_trios_35kbpos.txt", header=T)

# For each pedigree trio, evaluate whether the trio has inconsistent ancestry calls and if so, what is the inconsistent call
for (k in 1:nrow(trios)) {

        tmp <- trios[k,] # look at each pedigree trio separately
        
        # Get the column numbers in the ancestry calls data frame that correspond to the dad (d), mom (m), and kid (c) in the focal pedigree trio
        which(names(calls)==tmp$dad_id) -> d
        which(names(calls)==tmp$mom_id) -> m
        which(names(calls)==tmp$kid_id) -> c

        # Get chromosome, position, and ancestry call columns for all individuals in the focal pedigree trio
        tmp2 <- calls[c(1:2, d,m,c)]
        
        # Per position, count the number of individuals in the focal pedigree trio where that site fell on the exact base pair that transitioned between ancestry states (denoted by a "boundary" in 2get_ancestry_calls script) 
        tmp2$boundary <- rowSums(tmp2[c(3:5)]=="boundary")
        # Per position, count the number of individuals in the focal pedigree trio where the site did not fall in a tract because the tract was removed because it was <1 kb (dnoete dby a "not_in_tract" in 2get_ancestry_calls script)
        tmp2$not_in_tract <- rowSums(tmp2[c(3:5)]=="not_in_tract")
        # Per position, count the number of individuals in the focal pedigree trio that are heterozygous (i.e., ancestry state = 1 so one anubis allele and thus also one yellow allele)
        #tmp2$hetero_count <- rowSums(tmp2[c(3:5)]==1)
        
        # If all individuals in the focal pedigree trio have an ancestry state call (i.e., they do not have "boundary" or "not_in_tract" designation), define position as "usable", otherwise define as "not_usable" (since we need calls for all individuals in the trio to identify whether the ancestry calls are consistent in the pedigree trio)
        tmp2$usable_sites <- ifelse(tmp2$boundary==0 & tmp2$not_in_tract==0, "usable", "not_usable")

        # 1 is inconsistent, 0 is not inconsistent and trio is usable, NA is trio is not usable due to not having a call for one of the individuals (boundary, not in tract)
        tmp2$inconsistent <- ifelse(tmp2$usable_sites=="not_usable", NA, # assign NA to not usable sites for trios
        ifelse(tmp2$usable_sites=="usable" & as.character(tmp2[,3])==as.character(tmp2[,4]) & as.character(tmp2[,3])!="1" & as.character(tmp2[,5])!=as.character(tmp2[,3]), 1, # if parents are homozygous for the same ancestry and the offspring is not homozygous for the homozygous ancestry of its parents, assign site as inconsistent (1)
        ifelse(tmp2$usable_sites=="usable" & (as.character(tmp2[,3])==0 & as.character(tmp2[,4])==2 |  as.character(tmp2[,3])==2 & as.character(tmp2[,4])==0) & as.character(tmp2[,5])!=1, 1, # if parents are homozygous for opposite ancestries and the offspring is not heterozygous, assign site as inconsistent (1)
        ifelse(tmp2$usable_sites=="usable" & (as.character(tmp2[,3])==0 & as.character(tmp2[,4])==1 | as.character(tmp2[,3])==1 & as.character(tmp2[,4])==0) & as.character(tmp2[,5])==2, 1, # if one parent is homozygous for yellow ancestry and the other parent is heterozygous, if the offspring is homozygous anubis, assign site as inconsistent (1)
        ifelse(tmp2$usable_sites=="usable" & (as.character(tmp2[,3])==2 & as.character(tmp2[,4])==1 | as.character(tmp2[,3])==1 & as.character(tmp2[,4])==2) & as.character(tmp2[,5])==0, 1, 0))))) # if one parent is homozygous for anubis ancestry and the other parent is heterozygous, if the offspring is homozygous yellow, assign site as inconsistent (1) and assign all other other sites as consistent (0)

        # If the focal pedigree trio is inconsistent and it is a case where there are two ways for the focal pedigree trio being inconsistent (i.e., when both parents are homozygous), get what erroneous ancestry is being called
        tmp2$inconsistent_error <- ifelse(tmp2$inconsistent==1 & as.character(tmp2[,3])==as.character(tmp2[,4]), as.character(tmp2[,5]), # if parents are homozygous for the same ancestry and there is an inconsistency in the offspring, get offspring's ancestry
        ifelse(tmp2$inconsistent==1 & as.character(tmp2[,3])!=as.character(tmp2[,4]) & as.character(tmp2[,3])!=1 & as.character(tmp2[,4])!=1, as.character(tmp2[,5]), NA)) # if parents are homozygous for opposite ancestries and there is an inconsistency in the offspring, get the offspring's ancestry       
        
        # Lastly, create a column containing the three ancestry calls per position per pedigree trio as a record (dad ancestry call, mom ancestry call, kid ancestry call in that order)
        tmp2$calls <- paste(as.character(tmp2[,3]), as.character(tmp2[,4]), as.character(tmp2[,5], sep="_"))

        # Name the column containing inconsistencies as the ids in the pedigree trio and a "_1"
        colnames(tmp2)[9] <- paste(colnames(tmp2)[3], colnames(tmp2)[4], colnames(tmp2)[5], "1", sep="_")
        # Name the column containing the inconsistent calls if there are two ways of being inconsistent as the ids in the pedigree trio and a "_2"
        colnames(tmp2)[10] <- paste(colnames(tmp2)[3], colnames(tmp2)[4], colnames(tmp2)[5], "2", sep="_")
        # Name the column containing the three ancestry calls per position as the ids in the pedigree trio and a "_3"
        colnames(tmp2)[11] <- paste(colnames(tmp2)[3], colnames(tmp2)[4], colnames(tmp2)[5], "3", sep="_")      
       
        tmp3 <- tmp2[c(1:2,9)] # grab chromosome, position, and inconsistent columns
        tmp4 <- tmp2[c(1:2,10)] # grab chromosome, position, and inconsistent error columns
        tmp5 <- tmp2[c(1:2,11)] # grab chromosome, position, and calls columns

        if (k==1) {tmp3 -> inconsistent} # if this is the first pedigree trio, create a new data frame called "inconsistent"
        if (k>1) {merge(inconsistent, tmp3, by=c("chrom", "pos")) -> inconsistent} # if this is not the first pedigree trio, add current pedigree trio data to previous pedigree trio data

        if (k==1) {tmp4 -> inconsistent_error} # if this is the first pedigree trio, create a new data frame called "inconsistent_error"
        if (k>1) {merge(inconsistent_error, tmp4, by=c("chrom", "pos")) -> inconsistent_error} # if this is not the first pedigree trio, add current pedigree trio data to previous pedigree trio data
        
        if (k==1) {tmp5 -> calls_record} # if this is the first pedigree trio, create a new data frame called "calls_record"
        if (k>1) {merge(calls_record, tmp5, by=c("chrom", "pos")) -> calls_record} # if this is not the first pedigree trio, add current pedigree trio data to previous pedigree trio data


  print(k) # keep track of progress - print the pedigree trio number that has just finished       

}

# Save all three data frames into a single R Data file
save(inconsistent, inconsistent_error, calls_record, file="inconsistent_local_ancestry_calls_pedigree_trios_maskedSNPRCref.Rd")
print("done")
q(save="no")
