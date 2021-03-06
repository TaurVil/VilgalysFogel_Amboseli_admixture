# Panubis1Liftover
Generating a liftover file between _Panubis1_ and the human genome. These scripts are based on the UCSC pipelines.

# Sets of scripts
The process to create a liftover file is divided into 4 scripts, run in sequential order. The second script calls a subscript which is also included in this directory. These scripts were initially written for the UCSF cluster, which has some different systems than HARDAC at Duke. Some tweaking might be necessary to get everything to run smoothly on your specific computing cluser.

The files 1-4 refer to instructions to make the liftover from _hg38_ to _Panubis1_. The 1b-4b are files to make the liftover for _Panubis1_ to _hg38_. 
 
Note that both genomes are supposed to be repeat masked first, which they are when downloaded from NCBI (soft-masking is okay, which is when repetitive regions are in lower case). I believe all the necessary executables can be found here: http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/.
 
