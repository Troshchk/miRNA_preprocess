#!/bin/bash


PATH1="/home/ksenia/diploma/real/"
#PATH1="/home/ksenia/diploma/real"
FILES=$(ls $PATH1/*_first_fastqc/*.html | grep -o '..................$' | cut -c -6)

function multiPerSample() {
	PATH1="/home/ksenia/diploma/real"
	echo $1
	mkdir $PATH1/${1}_multi_per_sample
	cp $PATH1/${1}_first_fastqc/${1}_fastqc.html $PATH1/${1}_multi_per_sample/
	cp $PATH1/${1}_first_fastqc/${1}_fastqc.zip $PATH1/${1}_multi_per_sample/
	fastqc $PATH1/${1}_cutadapt/${1}_trimmed.fastq -o $PATH1/${1}_multi_per_sample/
	fastqc $PATH1/${1}_mapped/${1}*_mapped.fastq -o $PATH1/${1}_multi_per_sample/
	fastqc $PATH1/${1}_sortme/${1}_trimmed_norRNA_notRNA.fastq -o $PATH1/${1}_multi_per_sample/
	multiqc $PATH1/${1}_multi_per_sample/ -o $PATH1/${1}_multi_per_sample/
}

function multiPerSampleTest() {
	PATH1="/home/ksenia/diploma/real"
	echo $PATH1
	echo $1
	#mkdir $PATH1/$1_multi_per_sample
	#cp $PATH1/*_first_fastqc/*.html $PATH1/$1_multi_per_sample/
	echo $PATH1/$1_first_fastqc/$1_fastqc.html
	echo $PATH1/$1_cutadapt/$1_trimmed.fastq $PATH1/$1_multi_per_sample/
	echo $PATH1/${1}_mapped/${1}*_mapped.fastq  $PATH1/${1}_multi_per_sample/
	echo $PATH1/$1_sortme/${1}_trimmed_norRNA_notRNA.fastq  $PATH1/$1_multi_per_sample/
	echo $PATH1/$1_multi_per_sample/
}


#echo $PATH1
#for f in $FILES; do multiPerSampleTest $f; done 
export -f multiPerSampleTest
#parallel -j 16 multiPerSampleTest ::: $FILES
#mkdir $PATH1/1
export -f multiPerSample
#multiPerSample /47442S
parallel -j 16 multiPerSample ::: $FILES

