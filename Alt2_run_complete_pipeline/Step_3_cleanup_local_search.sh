#!/bin/bash
#SBATCH --partition=debug
#SBATCH -J step1_setupPA
#SBATCH -o setupPA.log
#SBATCH -t 1:00:00
#SBATCH -n 1
#SBATCH -c 1

mkdir shfiles
mkdir pickles
mkdir logs
mv *local1.pickle pickles/
mv *localSearch.sh shfiles/
mv *.log logs/

module load intelpython/2.7.12
make
python parallel_local_search_find_node_importance.py