#!/bin/bash
for i in *.fastq
do 
	i2=${i::-6}$"_trim.fastq"
	echo $i2
	java -jar /opt/Trimmomatic-0.35/trimmomatic.jar SE $i $i2 ILLUMINACLIP:/opt/Trimmomatic-0.35/adapters/TruSeq3-SE.fa:2:30:10 MINLEN:20
done
