#!/bin/bash

######################################################################################
### RUNNING THE FIRST FASTQC, CREATING FOLDER 'SAMPLE_FIRST_FASTQC' FOR THE OUTPUT ###
######################################################################################
	

function firstFastQC(){
	#echo 'RUNNING FIRST FASTQC' $1
	echo 'OUTPUT FOLDER' ${2}_first_fastqc
	mkdir ${2}_first_fastqc
	fastqc $1 -o ${2}_first_fastqc
}	



######################################################################################
### UNZIP AND EXTRACT THE SEQUENCES FOR TRIMMING: ADAPTORS AND OVERREPRESENTED #######
########################### WITH 'NO HIT' ############################################
######################################################################################


function extract(){

echo 'EXTRACTING SEQUENCES FOR CUTTING' $1
unzip ${2}_first_fastqc/*_fastqc.zip -d ${2}_first_fastqc #unzipping the fastqc results
mv ${2}_first_fastqc/${3}_fastqc ${2}_first_fastqc/$3 #renaming the output folder
#extracting all no hit seq from the overrepresented seq 
sed -n '/>>Overrepresented/,/>>END_MODULE/p' ${2}_first_fastqc/$3/fastqc_data.txt | grep -v 'No Hit\|>>\|#'  | awk '{print $1}' > ${2}_first_fastqc/$3/tocut.txt  
#extracting data for the adapter content from the fastqc report
sed -n '/>>Adapter/,/>>END_MODULE/p' ${2}_first_fastqc/$3/fastqc_data.txt | sed -n '/^[0-9]/p' >  ${2}_first_fastqc/$3/adapter.txt 
#cat $FOLDER/adapter.txt
#in case the adapter content is higher, than 0.1 and the number does not contain E (so not too small) the corresponding adapter identifiers are written to the adapter2.txt file
cat ${2}_first_fastqc/$3/adapter.txt | while read LINE; do  awk '{ if ($2 > 0.1 && ! /E/) {print 1}}' $line >> ${2}_first_fastqc/$3/adapter2.txt; done
cat ${2}_first_fastqc/$3/adapter.txt | while read LINE; do  awk '{ if ($3 > 0.1 && ! /E/) {print 2}}' $line >> ${2}_first_fastqc/$3/adapter2.txt; done
cat ${2}_first_fastqc/$3/adapter.txt | while read LINE; do  awk '{ if ($4 > 0.1 && ! /E/) {print 3}}' $line >> ${2}_first_fastqc/$3/adapter2.txt; done
cat ${2}_first_fastqc/$3/adapter.txt | while read LINE; do  awk '{ if ($5 > 0.1 && ! /E/) {print 4}}' $line >> ${2}_first_fastqc/$3/adapter2.txt; done
cat ${2}_first_fastqc/$3/adapter.txt | while read LINE; do  awk '{ if ($6 > 0.1 && ! /E/) {print 5}}' $line >> ${2}_first_fastqc/$3/adapter2.txt; done
#adapter2.txt sorted for the uniques adapter identifiers
sort ${2}_first_fastqc/$3/adapter2.txt | uniq > ${2}_first_fastqc/$3/adapter3.txt
#adapter sequences are added to the tocut.txt file for the further trimming
if grep -q '1' --max-count=1 ${2}_first_fastqc/$3/adapter3.txt; then echo "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT" >> ${2}_first_fastqc/$3/tocut.txt; fi
if grep -q '2' --max-count=1 ${2}_first_fastqc/$3/adapter3.txt; then echo "TGGAATTCTCGGGTGCCAAGG" >> ${2}_first_fastqc/$3/tocut.txt; fi
if grep -q '3' --max-count=1 ${2}_first_fastqc/$3/adapter3.txt; then echo 'GATCGTCGGACT' >> ${2}_first_fastqc/$3/tocut.txt; fi
if grep -q '4' --max-count=1 ${2}_first_fastqc/$3/adapter3.txt; then echo 'CTGTCTCTTATACACATCT' >> ${2}_first_fastqc/$3/tocut.txt; fi
if grep -q '5' --max-count=1 ${2}_first_fastqc/$3/adapter3.txt; then echo 'CGCCTTGGCCGT' >> ${2}_first_fastqc/$3/tocut.txt; fi
#all the temporary files are removed
rm ${2}_first_fastqc/$3/adapter*.txt
#each second line is appended as >adapter in order to get tocut.txt to the fasta format for cutadapt
sed -i '0~1s/^/>adapter\n/' ${2}_first_fastqc/$3/tocut.txt
}



######################################################################################
####### CUTADAPT FOR SOMAGENICS: trim adaptors/primers and cut the length at 15 ######
######################################################################################

function trim(){
echo 'CUTTING ADAPTORS' $1
ADAPTERS=${2}_first_fastqc/$3/tocut.txt
OUTPUT=${2}_cutadapt/${3}_trimmed.fastq.gz
UN_OUTPUT=${2}_cutadapt/${3}_without_adaptor.fastq.gz
LOG=${2}_cutadapt/${3}.log
#echo $ADAPTERS $OUTPUT $UN_OUTPUT $LOG
mkdir ${2}_cutadapt
#cutting adaptors from the given file containing adaptors and overrepresented sequences, the minimum length is 15, -n 5 for trimming all the adaptors each read is processed 5 times
cutadapt -a file:$ADAPTERS -o $OUTPUT --untrimmed-output $UN_OUTPUT -n 5 -m 35 -q 20 $FILE > $LOG
}


		
######################################################################################
################################ SORTME ##############################################
######################################################################################

#gunzip for decompressing, as the sortme2 is used and it does not work for the compressed files
#--fastx for fasta output
#--num-alignments 1: very fast, output the first alignment passing E-value, best choice for filtering only
#--other contains unaligned reads, so the ones, which are NOT rRNAs, resp. tRNAs -> what is needed
#rRNA references are build-in , for tRNA references:
#rsync -avzP rsync://ftp.ensembl.org/ensembl/pub/release-96/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.toplevel.fa.gz .
#rsync -avzP rsync://ftp.ensembl.org/ensembl/pub/release-96/gff3/homo_sapiens/Homo_sapiens.GRCh38.96.gff3.gz .
#rsync -avzP rsync://ftp.ensembl.org/ensembl/pub/release-96/gff3/homo_sapiens/Homo_sapiens.GRCh38.96.abinitio.gff3.gz .
#rsync -avzP rsync://ftp.ensembl.org/ensembl/pub/release-96/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz . 
#zcat Homo_sapiens.GRCh38.96.gff3.gz |grep -i trnascan > tRNAs.gff
#gunzip Homo_sapiens.GRCh38.dna_sm.toplevel.fa.gz
#bedtools getfasta -s -fi ./Homo_sapiens.GRCh38.dna_sm.toplevel.fa -bed ./tRNAs.gff -fo ./tRNAs.fasta
#indexdb_rna --ref /home/ksenia/diploma/sortme/tRNA/tRNAs.fasta,./tRNA/tRNAs: -v &> ./tRNA/indexdb_rna.log



function sortOut(){
echo 'SORTING' $1
#echo 'OUTPUT FOLDER' ${2}_sortme
#mkdir ${2}_sortme
#gunzip ${2}_cutadapt/${3}_trimmed.fastq.gz -d ${2}_cutadapt/
#/home/kseniatrs/external/sortmerna/bin/sortmerna --ref /home/kseniatrs/external/sortmerna/DBs/silva-euk-18s-id95.fasta --ref /home/kseniatrs/external/sortmerna/DBs/silva-euk-28s-id98.fasta --ref /home/kseniatrs/external/sortmerna/DBs/rfam-5s-database-id98.fasta --ref /home/kseniatrs/external/sortmerna/DBs/rfam-5.8s-database-id98.fasta --fastx --num_alignments 1 --reads $1 --aligned ${2}_sortme/${3}_trimmed_rRNA --other ${2}_sortme/${3}_trimmed_norRNA --workdir ${2}_sortme

mkdir ${2}_sortmet
/home/kseniatrs/external/sortmerna/bin/sortmerna --ref /home/kseniatrs/external/sortmerna/DBs/tRNAs.fasta --reads ${2}_sortme/${3}_trimmed_norRNA.fastq --aligned ${2}_sortmet/${3}_trimmed_norRNA_tRNA --fastx --other ${2}_sortmet/${3}_trimmed_norRNA_notRNA --workdir ${2}_sortmet
}





######################################################################################
################################ SALMON ##############################################
######################################################################################






######################################################################################
################################ RUNNING #############################################
######################################################################################


PATH1="/home/kseniatrs/external/2020_02_Prucha/splits" #IDENTIFY THE FOLDER WITH THE INPUT DATA
PATH2="/home/kseniatrs/external/2020_02_Prucha/splits"
FILES=$(find $PATH1 -maxdepth 1 -name '*.fastq.gz')

function finRun(){
	FILE=$1
	#echo 'FILE' $FILE
	FOLDER=${FILE%.fastq.gz}
	#echo 'OLD FOLDER' $FOLDER
	FOLDER=${FOLDER/$PATH1/$PATH2}
	BASE=${FOLDER##*/}
	DIR=$BASE
	echo 'NEW FOLDER' $FOLDER
	echo $BASE $DIR
	#firstFastQC "$FILE" "$FOLDER" "$BASE"
	#extract "$FILE" "$FOLDER" "$BASE"
    	#trim  "$FILE" "$FOLDER" "$BASE"
    	sortOut "$FILE" "$FOLDER" "$BASE"
    	#map "$FILE" "$FOLDER" "$BASE"
    	#echo 'DONE FOR ' $FILE
	#find -type f -name "$PATH2/*/*.fastq" -exec gzip "{}" 
	#for i in `ls $PATH2/*_mapped/*.sam`; do sample_name=`echo $i | sed "s/sam/bam/g"`; samtools view -Sb $i > $sample_name; done
}


#finRun '/home/ksenia/diploma/real/46634S.fastq.gz'
#for f in $FILES; do finRun $f; done
#export -f firstFastQC
#export -f extract
#export -f UMI
#export -f trim
export -f sortOut
#export -f map
export -f finRun
parallel  finRun ::: $FILES

