# This directory contains two scripts for evaluating whether genomic features predict the rate of change in anubis allele frequencies over time (Supplementary Methods 14.2).

#############################################################################################################################
# FST, mean recombination rate, and B values.
#############################################################################################################################

#############################################################################################################################
# Get pedigree inconsistences for the pedigree inconsistencies model (scripts labelled 1-3...).
#############################################################################################################################

# In R, run 1prep_for_pedigree_inconsistencies.R which generates two R data files (local_ancestry_pedigree_trios_maskedSNPRCref.Rd and local_ancestry_pedigree_trios_unmaskedWallref.Rd) containing tracts, pedigree trio info, list of pedigree individuals, and genomic positions for evaluating the consistent of ancestry calls within pedigree trios for two sets of ancestry tracts (one generated using the SNPRC reference panel, one generated using the Wall et al. 2016 Molecular Ecology low coverage reference panel). These R data files can then uploaded to a computing cluster for parallelization across chromosomes (which will make things run much faster).
