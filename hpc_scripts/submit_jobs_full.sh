path=/wehisan/home/allstaff/t/tonkin-hill.g/global_var_manuscript/mosaic_processed_data/
for FILE in $path/*_run*.fasta;
do
        qsub -F "$FILE" torque_run_full_iteration.sh
done
