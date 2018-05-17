path=/wehisan/home/allstaff/t/tonkin-hill.g/global_var_manuscript/mosaic_processed_data/
for NUM in 285 332 339 352 375 617 618;
do
        FILE=${path}Protein_NoLab_translateable_combined_454_tessema_centroids_reducedRegion_run${NUM}.fasta;
        qsub -F "$FILE" torque_run_3rd_viterbi_iteration.sh
done

