#!/bin/bash
chmod -R 755 MGP_pyGAsubmit.sh
for graphfilename in *.gpickle; do
	chmod -R 755 $graphfilename;
		for datfilename in *.bin; do
			chmod -R 755 $datfilename;
				sbatch MGP_pyGAsubmit.sh $graphfilename $datfilename;
		done		
done
