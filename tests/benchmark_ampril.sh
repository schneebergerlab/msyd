#!/bin/sh
# small script that 
cd /home/lrauschning/pansyri/pansyri
#git pull && python setup.py install
curtime=$(date --rfc-3339=seconds)
commit=$(git show --oneline -s) 
cores=4
cd /home/lrauschning/pansyri/data/ampril

rawout=''
for x in seq 1 10
do
	rawout=$rawout$(/usr/bin/time -f "%e,%M;" pansyri call -i full.tsv -v /dev/null -c $cores 2>&1 > /dev/null | tail -n 1)
done
echo $rawout
avgs=$(echo "rawout='$rawout'.split(';')[:-1];print('time:', ' s, mem: '.join([str(sum([float(tup.split(',')[ind]) for tup in rawout])/len(rawout)) for ind in [0, 1]]), 'Kb')" | python) # it's not ugly if it works!

# in case of not doing avgs:
#avgs=$(/usr/bin/time -f "time: %e s, mem: %M Kb" pansyri call -i full.tsv -v /dev/null -c $cores 2>&1 > /dev/null | tail -n 1)

echo $curtime – $commit: $avgs on $cores cores at $(hostname) >> /home/lrauschning/pansyri/pansyri/benchmarks.txt
