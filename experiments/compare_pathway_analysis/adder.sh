for i in `seq 0 5 20`; do
	touch exp_raw_${i}.txt
	for j in {0..9}; do
		sed -n 's/0,04350,/ /p' ${i}_${j}/pvalues2.csv >> exp_raw_${i}.txt
	done
done

for i in `seq 0 5 20`; do
	touch exp_adj_${i}.txt
	for j in {0..9}; do
		sed -n 's/0,04350,/ /p' ${i}_${j}/pvalues3.csv >> exp_adj_${i}.txt
	done
done


for i in `seq 0 5 20`; do
	touch log_raw_${i}.txt
	for j in {0..9}; do
		sed -n 's/0,04350,/ /p' ${i}_${j}/pvalues4.csv >> log_raw_${i}.txt
	done
done


for i in `seq 0 5 20`; do
	touch log_adj_${i}.txt
	for j in {0..9}; do
		sed -n 's/0,04350,/ /p' ${i}_${j}/pvalues5.csv >> log_adj_${i}.txt
	done
done


for i in `seq 0 5 20`; do
	touch mwu_${i}.txt
	for j in {0..9}; do
		sed -n 's/0,04350,/ /p' ${i}_${j}/pvalues6.csv >> mwu_${i}.txt
	done
done
