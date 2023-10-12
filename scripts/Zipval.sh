#!/bin/bash
module load anaconda/2.5

#file = basename of the file
#Input_File = Input Library
#Figure_Folder = where you store the generated figures

for file in ./*.gz
do	
	python calc_zipval.py $file Input_File Figure_Folder > $file.zipval.csv
gzip *.csv
done
