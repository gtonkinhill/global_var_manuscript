#!/bin/bash
#
#PBS -q u2296-1
#PBS -l nodes=1:ppn=1              # number of nodes
#PBS -o ./errors/slurm.$PBS_JOBID.out        # STDOUT
#PBS -e ./errors/slurm.$PBS_JOBID.err        # STDERR

cd /home/tonkin-hill.g/full_mosaic_run/results_iter3/

fasta=`basename $1`
/home/users/allstaff/tonkin-hill.g/global_var_manuscript/mosaic/mosaic -seq $1 -aa -tag $fasta -group 2 db target -target target -del 0.0114157102622 -eps 0.257764059982 -rec 0.0 > /home/tonkin-hill.g/full_mosaic_run/results_iter3/${fasta}_output.log

