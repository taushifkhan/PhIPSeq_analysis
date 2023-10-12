#!/bin/bash
module load anaconda/2.5
module load python/2.7

#file = your input file
#Metadata = virscan annotation file
#Beads = beadsnhit file
#Samples = samplesnhit file
#GROUPING_LEVEL = organism or species

for file in ./*.gz
do 
python calc_scores.py $file Metadata Beads Samples GROUPING_LEVEL 7 > $file.zscore.csv
gzip $file.zscore.csv
done
