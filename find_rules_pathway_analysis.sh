#!/bin/bash
for graphfilename in *.gpickle; do
	chmod -R 755 $graphfilename;
	sbatch calcNodeImportancesubmit.sh $graphfilename;
done