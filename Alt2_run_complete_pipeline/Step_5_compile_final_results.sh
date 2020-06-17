#!/bin/bash
#SBATCH --partition=debug
#SBATCH -J step1_setupPA
#SBATCH -o setupPA.log
#SBATCH -t 1:00:00
#SBATCH -n 1
#SBATCH -c 1

module load intelpython/2.7.12
make
python pathway_analysis_score_pathways_parallel_local_search.py  -sep , "dataName.csv" "matrix.name" "diff.file"