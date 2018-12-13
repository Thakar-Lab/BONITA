
for i in `seq 0 5 20`; do
	for j in {0..9}; do 
		cp matrix.csv ${i}_$j/
		cp diffFile.txt ${i}_$j/
		cp cleanup.sh ${i}_$j/
		cd ${i}_$j/
		bash cleanup.sh
		module load intelpython/2.7.12
		echo hsa04350_${i}_true_${j}.csv
		python pathway_analysis_score_pathways.py hsa04350_${i}_true_${j}.csv matrix.csv diffFile.txt -sep , > setup.txt
		cd ..
	done
done
