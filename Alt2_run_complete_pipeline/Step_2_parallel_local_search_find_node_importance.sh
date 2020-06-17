#!/bin/bash
#SBATCH --partition=debug
#SBATCH -J step1_setupPA
#SBATCH -o setupPA.log
#SBATCH -t 1:00:00
#SBATCH -n 1
#SBATCH -c 1

module load intelpython/2.7.12
make

for graphfilename in *.gpickle; do
	chmod -R 755 $graphfilename;
	for iteration in $(seq 1 5 ); do 
   		python re_run_nodes.py $graphfilename $iteration
   	done
done
