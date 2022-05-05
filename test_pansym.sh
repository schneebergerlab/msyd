#!/bin/sh
# full list of plants
# eri sha kyo col c24 ler an1 cvi
ref=an1
others="c24 col cvi eri kyo ler sha"
echo "usind as reference: $ref"
for i in $others
do
	for j in $others
	do
		#echo "$i & $j"
		#/usr/bin/time -f "%E" python ../pansr/pansym.pyx ampril/${ref}_${i}syri.out ampril/${ref}_${j}syri.out
		#continue
		for k in $others
		do
			echo "$i & $j & $k"
			/usr/bin/time -f "%E" python ../pansr/main.py ampril/${ref}_${i}syri.out ampril/${ref}_${j}syri.out ampril/${ref}_${k}syri.out
		done
	done
done

