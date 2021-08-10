#!/bin/bash
#SBATCH --get-user-env

g=${SLURM_ARRAY_TASK_ID}
j=`ls *.fa | head -$g | tail -1`

module load blat

## so we're using a rather lenient minimum similarity. We can edit that later, but I know that 98% is too stringent
blat ${MAINDIR}/fasta_old/${OLD%.fa*}.2bit ${j} ${j%.fa}_temp.psl -ooc=${OOC} -tileSize=${tileSize} -dots=10 -minScore=100 -minIdentity=90 -q=dna

