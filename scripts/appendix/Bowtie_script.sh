#!/bin/bash
module load Bowtie2/2.2.6
module load SAMtools/1.6

#file is the basename of the file
#genome = path to reference genome files
for file in ./*.fastq
	do 
		bowtie2 --very-sensitive-local -x genome $file -S $file.sam
		samtools view -Sb $file.sam > $file.bam
		samtools sort $file.bam > $file.sorted.bam
		samtools index $file.sorted.bam
		samtools idxstats $file.sorted.bam | cut -f 1,3 | sed -e '/^\*\t/d' -e "1 i id\file" | tr "\\t" "," > $file.count.csv
		gzip f$file.count.csv
done

