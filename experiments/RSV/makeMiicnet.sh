#!/bin/sh
#SBATCH --partition=standard
#SBATCH -J GArun
#SBATCH -o GA_out_%A_%a.txt
#SBATCH -t 36:00:00
#SBATCH -n 1

module load intelpython/2.7.12
srun python make_miic_net.py -t rsv.extreme.phenotype.CD4.tophat.rpm.filtered.20150310.txt connected_cd4_miic.graphml
echo "ran GA"