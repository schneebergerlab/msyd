#!/bin/sh
# full list of plants
# eri sha kyo col c24 ler an1 cvi
ref=cvi
echo "usind as reference: $ref"
for i in eri sha kyo col c24 ler an1;
do
	for j in eri sha kyo col c24 ler an1
	do
		echo "$i on $j"
		/usr/bin/time -f "%E" python ../pansr/pansym.pyx ampril/${ref}_${i}syri.out ampril/${ref}_${j}syri.out
	done
done

