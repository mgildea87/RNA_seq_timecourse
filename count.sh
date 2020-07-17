#!/bin/bash
A_filepath="/home/mag456/BioMG8310_2017/Fastq_files/"
D_file_path="/home/mag456/BioMG8310_2017/aligned_counted_concat/"
cere_gff="/home/mag456/BioMG8310_2017/genome/cere_gen.gff"
pom_gff="/home/mag456/BioMG8310_2017/genome/sp_feature.gff3"
con_cat_gff="/home/mag456/BioMG8310_2017/genome/concat_features.gff"
for i in *.bam
do
	sample=${i::-4}
	echo $sample 
	echo $i
	python -m HTSeq.scripts.count -f bam -s reverse -t gene -m intersection-nonempty -i ID $i $con_cat_gff > $D_file_path$sample"_concat_counts.txt"
	#python -m HTSeq.scripts.count -f bam -t gene -m intersection-nonempty $i $pom_gff > $D_file_path$sample"_pom_counts.txt"
done
