#!/bin/sh
#SBATCH --partition=preempt
#SBATCH -a 1-25
#SBATCH -J GArun
#SBATCH -o GA_out_%A_%a.txt
#SBATCH -t 24:00:00
#SBATCH -n 1

module load intelpython/2.7.12
srun python experiments.py $1 $SLURM_ARRAY_TASK_ID
echo "ran GA"