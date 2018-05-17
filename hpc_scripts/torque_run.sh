#!/bin/bash
#
#PBS -q u2296-4
#PBS -l nodes=1:ppn=1              # number of nodes
#PBS -o slurm.$PBS_JOBID.out        # STDOUT
#PBS -e slurm.$PBS_JOBID.err        # STDERR

cd /home/tonkin-hill.g/full_mosaic_run/results/

fasta=`basename $1`
/home/users/allstaff/tonkin-hill.g/mosaic_mod/mosaic -ma -seq $1 -aa -tag $fasta  -group 2 db target -target target -del 0.010761 -eps 0.396444 -rec 0.014 > /home/tonkin-hill.g/full_mosaic_run/results/${fasta}_output.log

