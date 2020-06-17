#!/bin/sh
#SBATCH --partition=standard
#SBATCH -a 1-10
#SBATCH -J ISrun
#SBATCH -o GA_out_%A_%a.txt
#SBATCH -t 48:00:00
#SBATCH -n 1
#SBATCH -c 1

module load intelpython/2.7.12
python calc_IS_parallel_local_search.py $1 $2 $SLURM_ARRAY_TASK_ID