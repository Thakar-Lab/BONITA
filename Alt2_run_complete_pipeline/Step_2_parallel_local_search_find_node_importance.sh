#!/bin/bash
#SBATCH --partition=standard
#SBATCH -a 1
#SBATCH -J GArun
#SBATCH -o Step2_node_importance_out_%a.txt
#SBATCH -t 48:00:00
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem 50G

module load intelpython/2.7.12
make

for graphfilename in *.gpickle; do
	chmod -R 755 $graphfilename;
    python parallel_local_search_find_node_importance.py $graphfilename $SLURM_ARRAY_TASK_ID
done
