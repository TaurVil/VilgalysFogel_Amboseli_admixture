# Panubis1Liftover
Generating a liftover file between Panubis1 and the human genome. These scripts are based on the UCSC pipelines

# Sets of scripts
There are two sets of scripts in this file. First, there is a large script "crossSpeciesLiftover" which contains the entire process. However, this is unwieldly to run, so we made and use the 4 "liftover_XX.sh" scripts for most practical computing. These scripts were initially written for the UCSF cluster, which has some different systems than HARDAC at Duke. Some tweaking might be necessary to get everything to run smoothly.

The files 1-4 refer to instructions to make the liftover from hg38 to panubis1. The 1b-4b are files to make the liftover for panubis1 to hg38. 
 
Note that both genomes are supposed to be repeat masked first, which they are when downloaded from NCBI (soft-masking is okay, which is when repetitive regions are in low case). I think all the necessary executables can be found here: http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/.
 
# Citation/Credit
Scripts were initially developed by Jacqueline Robinson based on UCSC procedures. They were later modified by Tauras Vilgalys and Jordan Anderson. Please contact me (Tauras) for access to the resulting recombination maps, with any questions, or for how to acknowledge/cite this code if used in your project. My e-mail is taur.vil@gmail.com. 
