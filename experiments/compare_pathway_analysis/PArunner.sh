mkdir trues
module load intelpython/2.7.12
for i in `seq 0 5 20`; do
	for j in {0..9}; do 
		mkdir ${i}_${j}
		cp *.py ${i}_$j/
		cp testRun.so ${i}_$j/testRun.so
		cp -r inputData ${i}_$j/inputData
		cp hsa04350.gpickle ${i}_${j}/hsa04350.gpickle
		cp hsa04350_${i}_true_${j}.csv ${i}_${j}/hsa04350_${i}_true_${j}.csv
		mv hsa04350_${i}_true_${j}.csv trues/hsa04350_${i}_true_${j}.csv
		cp calcNodeImportancesubmit.sh ${i}_$j/calcNodeImportancesubmit.sh
		cp just_tgf.gmt ${i}_$j/	
		cd  ${i}_${j}/
		python pathway_analysis_setup.py hsa04350_${i}_true_${j}.csv just_tgf.gmt -sep , > setup.txt
		for graphfilename in *.gpickle; do
			chmod -R 755 $graphfilename;
			sbatch calcNodeImportancesubmit.sh $graphfilename;
		done
		cd ..
	done
done 
