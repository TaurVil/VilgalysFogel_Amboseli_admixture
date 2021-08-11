## Ancestry and recombination rates 
library(data.table); library(ggplot2); library(parallel)

#########   Chunk for metadata       ###########
#### block size
d_name <- '5Mb'
distance <- 5e6

#########   Chunk to get data       ###########
### Recent vs historical
load("./windows.25kb.RData")

#########   Convert 25kb tracts to new tract lengths tracts       ###########
lengths <- read.delim("../03_baboon-genomic-resources/Panubis1_chromlengths.txt", header=F)
new_windows <- NULL; for (i in 1:20) {
  max <- lengths$V2[i]
  rows <- ceiling(max/distance)#; print(i); print(rows)
  tmp <- as.data.frame(matrix(ncol=3, nrow=rows))
  colnames(tmp) <- c("chr", "start", "end")
  tmp$chr <- paste("chr",i,sep="")
  tmp$start <- (0:(rows-1))*distance
  tmp$end <- (1:rows)*distance
  tmp <- subset(tmp, tmp$start >= 5e4 & tmp$end <= (max-5e4))
  rbind(new_windows,tmp) -> new_windows; rm(tmp, max, rows)
}; rm(i, lengths)

fun_get_new_window_info <- function(site, input_features) {
  tmp_features <- input_features[input_features$chr == new_windows$chr[site] & 
                                   input_features$start >= new_windows$start[site] & 
                                   input_features$end <= new_windows$end[site], ]
  
  return(c(colMeans(tmp_features[,4:5]), colSums(tmp_features[,6:17]), colMeans(tmp_features[,18:20])))
}

system.time(r <- mclapply(1:nrow(new_windows), FUN=fun_get_new_window_info, input_features=features))
r <- do.call("rbind",r)
cbind(new_windows,r) -> new_features

rm(r, fun_get_new_window_info, features, new_windows)

#### Filter to remove regions with recombination rate > 100x the median #####
max_rcr <- 100*median(new_features$recombination)
to_analyze <- subset(new_features, new_features$recombination < max_rcr)

##### Save product #####
save.image(paste("./windows.",d_name,".RData", sep=""))


##### Optional: get individaul ancestry per window ########
## this code scales `ancestry_per_individual`, if retained, to window sizes that match those above
# fun_get_new_window_ancestry <- function(site, input_ancestry) {
#   t_ancestry <- subset(input_ancestry, input_ancestry$chr == new_windows$chr[site] & input_ancestry$start >= new_windows$start[site] & input_ancestry$end <= new_windows$end[site])
#   
#   return(colMeans(t_ancestry[,4:(n+3)], na.rm=T)/2)
# }
# system.time(r <- mclapply(1:nrow(new_windows), FUN=fun_get_new_window_ancestry, input_ancestry=ancestry_full_1kb))
# window_by_individual <- as.data.frame(do.call("rbind",r)); rm(r)
# rm(new_windows, to_analyze, max_rcr, fun_get_new_window_info, distance_ancestry, distance_rcr, distance_B, fun_get_new_window_ancestry, rcr, B_25kb, fst_25kb, ancestry_full_1kb)
# save.image(paste("./",d_name,"_by_individual.RData", sep=""))

