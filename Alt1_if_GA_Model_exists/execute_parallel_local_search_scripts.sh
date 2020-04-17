#!/bin/sh
#SBATCH --partition=debug
#SBATCH -J localSearch
#SBATCH -o localSearch.log
#SBATCH -t 1:00:00
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=1G

for scriptname in *_localSearch.sh; do
	chmod -R 755 $scriptname;
	sbatch $scriptname;
done