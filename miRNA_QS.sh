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
############## UMI TOOLS FOR QIAGEN: extract UMIs and cut the length at 15 ###########
######################################################################################

function UMI(){
echo 'EXTRACTING UMIS' $1
LOG1=${2}_cutadapt/${3}.log
OUTPUT1=${2}_cutadapt/${3}_trimmed_long.fastq.gz
OUTPUT2=${2}_cutadapt/${3}_trimmed.fastq.gz
LOG2=${2}_cutadapt/${3}_trimmed_final.log
#echo $LOG1 $OUTPUT1 $LOG2 $OUTPUT2
mkdir ${2}_cutadapt	
umi_tools -v
#extracting UMIs with the pattern, which includes the adapter sequence (was found in documentation for UMI)
umi_tools extract --extract-method=regex --bc-pattern=".*(?P<discard_1>AACTGTAGGCACCATCAAT){s<=1}(?P<umi_1>.{10})(?P<discard_2>.*)" -L $LOG1 --stdin=$1  --stdout $OUTPUT1
#trimmimg length to min 15 bp (too small, could be broken miRNAs, rRNAs, another smallRNAs)
cutadapt -m 15 -q 20 -o $OUTPUT2 $OUTPUT1 > $LOG2
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
cutadapt -a file:$ADAPTERS -o $OUTPUT --untrimmed-output $UN_OUTPUT -n 5 -m 15 -q 20 $FILE > $LOG
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
mkdir ${2}_sortme
gunzip ${2}_cutadapt/${3}_trimmed.fastq.gz -d ${2}_cutadapt/
/home/ksenia/diploma/3.sortme/sortmeRNA/sortmerna --ref /home/ksenia/diploma/3.sortme/sortmeRNA/rRNA_databases/silva-euk-18s-id95.fasta,/home/ksenia/diploma/3.sortme/Indices/silva-euk-18s:/home/ksenia/diploma/3.sortme/sortmeRNA/rRNA_databases/silva-euk-28s-id98.fasta,/home/ksenia/diploma/3.sortme/Indices/silva-euk-28s:/home/ksenia/diploma/3.sortme/sortmeRNA/rRNA_databases/rfam-5s-database-id98.fasta,/home/ksenia/diploma/3.sortme/Indices/rfam-5s:/home/ksenia/diploma/3.sortme/sortmeRNA/rRNA_databases/rfam-5.8s-database-id98.fasta,/home/ksenia/diploma/3.sortme/Indices/rfam-5.8s --log --fastx --num_alignments 1 --reads ${2}_cutadapt/${3}_trimmed.fastq --aligned ${2}_sortme/${3}_trimmed_rRNA --other ${2}_sortme/${3}_trimmed_norRNA
/home/ksenia/diploma/3.sortme/sortmeRNA/sortmerna --ref /home/ksenia/diploma/3.sortme/tRNA/tRNAs.fasta,/home/ksenia/diploma/3.sortme/tRNA/tRNAs: --reads ${2}_sortme/${3}_trimmed_norRNA.fastq --log --aligned ${2}_sortme/${3}_trimmed_norRNA_tRNA --fastx --other ${2}_sortme/${3}_trimmed_norRNA_notRNA
}





######################################################################################
################################ BOWTIE ##############################################
######################################################################################

#mapping to the: 1.human mature miRNA, 2.human hairpin, 3.non-human mature, 4.non-human hairpin (non human might help to identify some conserved miRNAs) 
#wget -nc ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.gz
#wget -nc ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz
#divide to non/human sequences files hairpin and mature, build indexes
#bowtie-build mature_homo.fa mature_h > matureh.log
#bowtie-build hairpin_homo.fa hairpin_h > hairpinh.log
#bowtie-build mature_non_homo.fa mature_non_homo > mature_non_homo.log
#bowtie-build hairpin_non_homo.fa hairpin_non_homo > hairpin_non_homo.log


function map(){
echo 'MAPPING' $1
#echo 'OUTPUT FOLDER' ${2}_mapped
mkdir ${2}_mapped
bowtie --best --al ${2}_mapped/${3}_matureh_mapped.fastq --un ${2}_mapped/${3}_matureh_unmapped.fastq -t /home/ksenia/diploma/4.mapping/matureh_index/mature_h ${2}_sortme/${3}_trimmed_norRNA_notRNA.fastq -S ${2}_mapped/${3}_matureh_bowtie.sam --verbose &> ${2}_mapped/${3}_matureh_bowtie.log
bowtie --best --al ${2}_mapped/${3}_hairpinh_mapped.fastq --un ${2}_mapped/${3}_hairpinh_unmapped.fastq -t /home/ksenia/diploma/4.mapping/hairpinh_index/hairpin_h ${2}_mapped/${3}_matureh_unmapped.fastq -S ${2}_mapped/${3}_hairpinh_bowtie.sam --verbose &> ${2}_mapped/${3}_hairpinh_bowtie.log
bowtie --best --al ${2}_mapped/${3}_mature_non_h_mapped.fastq --un ${2}_mapped/${3}_mature_non_h_unmapped.fastq -t /home/ksenia/diploma/4.mapping/mature_non_homo/mature_non_homo ${2}_mapped/${3}_hairpinh_unmapped.fastq -S ${2}_mapped/${3}_mature_non_h_bowtie.sam --verbose &> ${2}_mapped/${3}_mature_non_h_bowtie.log
bowtie --best --al ${2}_mapped/${3}_hairpin_non_h_mapped.fastq --un ${2}_mapped/${3}_hairpin_non_h_unmapped.fastq -t /home/ksenia/diploma/4.mapping/hairpin_non_homo/hairpin_non_homo ${2}_mapped/${3}_mature_non_h_unmapped.fastq -S ${2}_mapped/${3}_hairpin_non_h_bowtie.sam --verbose &> ${2}_mapped/${3}_hairpin_non_h_bowtie.log
}


######################################################################################
################################ RUNNING #############################################
######################################################################################


PATH1="/home/ksenia/diploma/second_batch" #IDENTIFY THE FOLDER WITH THE INPUT DATA
PATH2="/home/ksenia/diploma/second_batch"
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
	firstFastQC "$FILE" "$FOLDER" "$BASE"
	extract "$FILE" "$FOLDER" "$BASE"
if [[ "$FILE" == *Q* ]]; then
	UMI  "$FILE" "$FOLDER" "$BASE"
else
    trim  "$FILE" "$FOLDER" "$BASE"
fi
    sortOut "$FILE" "$FOLDER" "$BASE"
    map "$FILE" "$FOLDER" "$BASE"
    #echo 'DONE FOR ' $FILE
	find -type f -name "$PATH2/*/*.fastq" -exec gzip "{}" 
	for i in `ls $PATH2/*_mapped/*.sam`; do sample_name=`echo $i | sed "s/sam/bam/g"`; samtools view -Sb $i > $sample_name; done
}


#finRun '/home/ksenia/diploma/real/46634S.fastq.gz'
export -f firstFastQC
export -f extract
export -f UMI
export -f trim
export -f sortOut
export -f map
export -f finRun
parallel -j 16 finRun ::: $FILES

