path=/wehisan/home/allstaff/t/tonkin-hill.g/global_var_manuscript/mosaic_processed_data/
for FILE in $path/*_run*.fasta;
do
fasta=`basename $FILE`
if [ -f "/home/tonkin-hill.g/full_mosaic_run/results_jump/${fasta}_align.txt" ]
then
  echo "$fasta found."
else
  qsub -F "$FILE" torque_run_full_iteration.sh
fi
done

