#!/bin/bash
for graphfilename in *.gpickle; do
	chmod -R 755 $graphfilename;
	sbatch pyGAsubmit.sh $graphfilename 1;
	sbatch pyGAsubmit.sh $graphfilename 2;
	sbatch pyGAsubmit.sh $graphfilename 3;
	sbatch pyGAsubmit.sh $graphfilename 4;
	sbatch pyGAsubmit.sh $graphfilename 5;
	sbatch pyGAsubmit.sh $graphfilename 6;
	sbatch pyGAsubmit.sh $graphfilename 7;
done