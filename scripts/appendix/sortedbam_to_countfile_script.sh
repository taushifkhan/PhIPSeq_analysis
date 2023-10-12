#!/bin/bash
module load Bowtie2/2.2.6
module load SAMtools/1.6

#location where the data is stored
data="/gpfs/projects/SDR_PPM1-1220-150017_nmarr/analysis/SidraOccHealth_vir3/run41_NS500_HKYTMBGX5/fastq_2018-XX-XX/PIDRUN3NGS48/vir3/data"
#new location where the meregd files will be saved
merged_data="/gpfs/projects/SDR_PPM1-1220-150017_nmarr/analysis/SidraOccHealth_vir3/run41_NS500_HKYTMBGX5/fastq_2018-XX-XX/PIDRUN3NGS48/pep2/merged"
#Change the location
genomeIndex="/gpfs/projects/SDR_PPM1-1220-150017_nmarr/bwh/reference_seq"
alignment="/gpfs/projects/SDR_PPM1-1220-150017_nmarr/analysis/SidraOccHealth_vir3/run41_NS500_HKYTMBGX5/fastq_2018-XX-XX/PIDRUN3NGS48/vir3/sam"
bamfiles="/gpfs/projects/SDR_PPM1-1220-150017_nmarr/analysis/SidraOccHealth_vir3/run41_NS500_HKYTMBGX5/fastq_2018-XX-XX/PIDRUN3NGS48/vir3/bam"
bamsorted="/gpfs/projects/SDR_PPM1-1220-150017_nmarr/analysis/SidraOccHealth_vir3/run41_NS500_HKYTMBGX5/fastq_2018-XX-XX/PIDRUN3NGS48/vir3/bamsorted"
countfiles="/gpfs/projects/SDR_PPM1-1220-150017_nmarr/analysis/SidraOccHealth_vir3/run41_NS500_HKYTMBGX5/fastq_2018-XX-XX/PIDRUN3NGS48/vir3/countfiles"
ARRAY=()

#reads all the files with extension fastq. Chnage the extension if the files are compressed to gz as *.gz or as *fastq.gz
for file in $data/*.fastq
	do      
        	#echo $file
                filename=${file##*/}
#               echo "filename:"$filename
##              basefilename=$(echo $filename | cut -d'.' -f -1)
                #echo $basefilename
#next line is based on what you had provided me as example files.
#SIR000037000788-IDX093-CCAAGTC_S9_L001_R1_001.fastq 
#filename will be split based on delimiter (-d'_') underscore and last two parts (-f -2) will be discarded.
#Result will be : SIR000037000788-IDX093-CCAAGTC_S9
#This is unique name for each sample so all the lanes (2 lanes or 4 lanes or 6 lanes,...) will be merged.  
                samplename=$(echo $filename | cut -d'_' -f -2)
#               echo $samplename

#There is a better way to do this but this should be ok for small list of files. Add each samplename to an ARRAY.
#for item in "${ARRAY[@]}"
#	do
#		[ "$samplename" != "$item" ] 
		ARRAY+=($samplename)
#echo $item

##              sampleid=$(echo $basefilename | cut -d'_' -f 2)
#               echo $sampleid
#done
done
#echo ${ARRAY[@]}
#After the for loop, sort the array and get the unique samplename
unique_ids=($(echo "${ARRAY[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))
#echo ${unique_ids[@]}

#merge lanes based on unique samplename and save it to a new folder and align as the files are merged.
for samplename in "${unique_ids[@]}"
do
#	cat $data/$samplename* > $merged_data/$samplename.fastq
#	echo "bowtie2 --very-sensitive -x $genomeIndex/vir2 $merged_data/$samplename.fastq -p 8 -S $alignment/${samplename}.sam"| bsub -J ${samplename} -P projectname -e err.txt -o out.log -n 8
# -w is the dependency argument in bsub. 
#	echo "samtools view -Sb $alignment/${samplename}.sam > $bamfiles/${samplename}.bam" | bsub -J ${samplename}_view -P projectname -e err.txt -o out.log -w "done('"${samplename}"')"
#change the -w to what you have named the previous -J such as -w "done('"${samplename}_view"')" and so .
#samtools sort
#   echo "samtools sort $bamfiles/${samplename}.bam -o $bamsorted/${samplename}.sorted.bam" | bsub -J ${samplename}_sort -P projectname -e err.txt -o out.log -w "done('"${samplename}_view"')"
#samtools index
#    echo "samtools index $bamsorted/${samplename}.sorted.bam" | bsub -J ${samplename}_index -P projectname -e err.txt -o out.log -w "done('"${samplename}_sort"')"
#samtools idx
echo "samtools idxstats $bamsorted/${samplename}_A.sorted.bam | cut -f 1,3 | sed -e '/^\*\t/d' -e '1 i id\t'"${samplename}"'' | tr '\t' , > $countfiles/${samplename}_count.csv" | bsub -J ${samplename}_idx -P virscan -e err.txt -o out.log



done

