#!/bin/bash
module load anaconda/2.5
module load python/2.7

#Figure_Folder: where you store the generated figures
#Threshold = significant p-val

for file in ./*.gz
do
	python call_hits.py $file-R1 $file-R2 Figure_Folder Threshold > file.zihit.csv
gzip *.csv
done
