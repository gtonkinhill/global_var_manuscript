#!/bin/bash
#
#PBS -q u2296-1
#PBS -l nodes=1:ppn=1              # number of nodes
#PBS -o ./errors/slurm.$PBS_JOBID.out        # STDOUT
#PBS -e ./errors/slurm.$PBS_JOBID.err        # STDERR

cd /home/tonkin-hill.g/full_mosaic_run/results_jump/

fasta=/wehisan/home/allstaff/t/tonkin-hill.g/global_var_manuscript/mosaic_processed_data/Protein_NoLab_translateable_combined_454_tessema_centroids_reducedRegion_randomSampleTargetSize1000.fasta

tag=Protein_NoLab_translateable_combined_454_tessema_centroids_reducedRegion_randomSampleTargetSize1000_jump_restart${1}

/home/users/allstaff/tonkin-hill.g/global_var_manuscript/mosaic/mosaic -seq $fasta -aa -tag $tag -group 2 db target -target target -del 0.0166765773591 -eps 0.273305573191 -rec $1 > /home/tonkin-hill.g/full_mosaic_run/results_jump/${tag}_output.log

