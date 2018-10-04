#!/bin/sh
#SBATCH --partition=standard
#SBATCH -a 1-10
#SBATCH -J GArun
#SBATCH -o GA_out_%A_%a.txt
#SBATCH -t 24:00:00
#SBATCH -c 1

module load intelpython/2.7.12
srun python experiments.py $1 $2 $SLURM_ARRAY_TASK_ID
echo "ran GA"