path=/wehisan/home/allstaff/t/tonkin-hill.g//mosaic_mod/full_run/combined_454_contamFree_translated0stops_NoLabIso_centroids_split/
for FILE in $path/*_run*.fasta;
do
	qsub -F "$FILE" torque_run.sh
done

