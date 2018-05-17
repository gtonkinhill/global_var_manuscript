path=/wehisan/home/allstaff/t/tonkin-hill.g/global_var_manuscript/mosaic_processed_data/
for jump in 0.001 0.013; 
#for jump in $(seq 0.00 0.001 0.1);
do
        qsub -F "$jump" torque_run_jump_viterbi_iterationi_restart.sh
done

