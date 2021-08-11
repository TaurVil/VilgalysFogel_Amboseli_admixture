#!/usr/bin/env Rscript
#SBATCH --get-user-env

library(data.table); library(plyr)
fread("NAME.output.txt", fill=T) -> data
read.delim("~/genomes/panubis1/Panubis1.0.fa.fai",header=F) -> lengths

# Column 1: number of generations for the output
# Column 2: 0, sub-population 0
# Column 3: 1/male; 0/female
# Column 4: 1 individual index, number of individuals printed
# Column 5: chromosome, 0 to n (8)
# Column 6: 0/1 (0 is maternal; 1 is paternal)
# Column 7: 0/1, ancestry type
# Column 8/9: position in M

output <- NULL
for (i in unique(data$V5)) {
        print(paste("starting chrom", i))

        d <- subset(data, data$V5 == i)

        # Let's clean up the file for that chromosome.
        d <- d[order(d$V9,d$V8),]
        # Remove tracks of no length (product of the simulation)
        d <- subset(d, !(d$V8 == d$V9))
        # Separate out maternal and paternal copies. Merge tracts that were joined by removing tracts of no length.
        d0 <- subset(d, d$V6 == 0); d1 <- subset(d, d$V6 == 1)
        d0$nxt <- c(d0$V7[-1],4);       d1$nxt <- c(d1$V7[-1],4)
        d0$nxt_pos <- c(d0$V9[-1],4); d1$nxt_pos <- c(d1$V9[-1],4)
        # Replace with the appropriate end point
        d0$V9[d0$V7 == d0$nxt] <- d0$nxt_pos[d0$V7 == d0$nxt]
        d1$V9[d1$V7 == d1$nxt] <- d1$nxt_pos[d1$V7 == d1$nxt]
        # Get ride of the second one in the string
        d0$delete <- c(0, (d0$V7 == d0$nxt)[-nrow(d0)]); d1$delete <- c(0, (d1$V7 == d1$nxt)[-nrow(d1)])
        d0 <- subset(d0, d0$delete == 0); d1 <- subset(d1, d1$delete == 0)
        # Remerge data from both chroms
        rbind(d0,d1) -> d; d <- d[order(d$V9,d$V8),]; rm(d1, d0)

        #Create and structure output file
        res <- as.data.frame(matrix(ncol=2,nrow=nrow(d)))
        res$chrom <- lengths[i+17,1]; res[,2] <- d$V9; res[,1] <- c(0,res[-nrow(res),2]); res$ancestry  <- NA
        for (j in 1:nrow(d)) {
                tmp  = subset(d, d$V9 == res$V2[j] |  (d$V9 >= res$V2[j]  & d$V8 < res$V2[j]) ); tmp
                # the tract that ends there | (or) the end is later
                nrow(tmp) -> res$rows[j]
                sum(tmp$V7) -> res$ancestry[j]
                #if (j %in% seq(from=10,to=1e7,by=100)) {print(j)}
        }
        res<-res[-nrow(res),]
        # Remove cases with no length (I don't think there will be any after our earlier changes)
        res<-subset(res, !(res$V1 == res$V2))
        # Remaining cases are where ancestry tracts switch between parental chromosomes. We'll fix those so this resembles our ancestry calls, but we should also note that it's rather rare to share the same break point.
        res$a_next <- c(res$ancestry[-1],-9) #; print(paste("repetitive calls ", sum(res$ancestry == res$a_next, na.rm=T),sep="_"))
        res$nxt_pos <- c(res$V2[-1],-9)
        res$V2[res$ancestry == res$a_next] <- res$nxt_pos[res$ancestry == res$a_next]

        res$delete <- c(0, (res$ancestry == res$a_next)[-nrow(res)])
        res <- subset(res, res$delete == 0)
        res$a_next <- c(res$ancestry[-1],-9); print(paste("repetitive calls ", sum(res$ancestry == res$a_next, na.rm=T),sep="_"))
        print(paste("more than 2 chroms:", sum(!res$rows == 2)))
        l <- lengths[i+17,2]
        m <- max(res$V2)
        res$V1 <- round(res$V1*l/m )
        res$V2 <- round(res$V2*l/m )
        print(paste("min_length", min(res$V2 - res$V1)))

        rbind(output, res) -> output ; rm(res,j, tmp, d )
}; rm(i)
colnames(output)[1:2] <- c('start','end')
write.table(output[,1:4], "tracts/NAME.tracts.txt",row.names=F, col.names=T, sep="\t", quote=F)

