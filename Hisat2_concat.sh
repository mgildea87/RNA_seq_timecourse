#!/bin/bash
for i in *_trim.fastq
do 
        outname=${i::-11}
	echo $outname
	hisat2 -p 4 --phred33 --no-unal -x ~/BioMG8310_2017/genome/concat_genome_index --max-intronlen 1000 -U $i | samtools view -bhq 10 - | samtools sort - -o $outname.bam 
done
