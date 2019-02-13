#!/bin/bash

for dataset in *_data; do
	chmod -R 755 $dataset;
	cd $dataset
	for dataset in *_data.csv; do	
		module load intelpython/2.7.12
		python pathway_analysis_setup.py $dataset True_positive_diseases.gmt > setup.txt
	done	
	for graphfilename in *.gpickle; do
		chmod -R 755 $graphfilename;
		sbatch calcNodeImportancesubmit.sh $graphfilename;
	done
	cd ..
done
