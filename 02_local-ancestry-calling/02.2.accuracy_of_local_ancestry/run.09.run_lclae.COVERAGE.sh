#!/bin/bash

index=${SLURM_ARRAY_TASK_ID}

id=$index
cov=COVERAGE
for chrom in `echo CHROMOSOME`; 
	do name=`head -$id tmp.012.indv | tail -1 | sed -e s/.varia.*$//`; 
		for d in `cat 02_AIMs.txt`; do for window in `cat 02_windows.txt`; 
				do touch ./raw_calls/full_refpanel.${window}.${d}.$name.$cov.$chrom.txt; 
				
				l=`wc -l raw_calls/full_refpanel.${window}.${d}.$name.$cov.$chrom.txt | cut -f 2 | sed 's/ .*//g'`; 
				
				if (($l==0)); then echo "STARTING"; 
				
				grep ^$chrom genolik.$cov.genolik | ~/Programs/LCLAE/filtbaboon2c 56 anubis.h yellow.h $id | ~/Programs/LCLAE/geno_lik2 ${d} ${window} > ./raw_calls/${window}.${d}.$name.$cov.$chrom.txt; fi; 
				
echo $id $cov $chrom $d $window; done; done; done

