#!/bin/bash
#SBATCH --partition=debug
#SBATCH -J step1_setupPA
#SBATCH -o setupPA.log
#SBATCH -t 1:00:00
#SBATCH -n 1
#SBATCH -c 1

module load intelpython/2.7.12
python pathway_analysis_setup_net.py  -sep , "dataName.csv" "pathway.graphml"

for graphfilename in *.gpickle; do
	chmod -R 755 $graphfilename;
	sbatch calcNodeImportancesubmit_parallel_local_search.sh $graphfilename;
done