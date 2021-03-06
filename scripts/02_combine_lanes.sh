#!/bin/bash 
#SBATCH --job-name=combine_lanes
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=5G
#SBATCH --qos=general
#SBATCH --partition=general
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err



# combine data across lanes for each sample/mate pair

# input directories
SMALL=../data/small_fastq
TOTAL=../data/total_fastq

# create output directories
COMBS=../data/small_combined
if [ -d $COMBS ]; then rm -r $COMBS; fi
mkdir -p $COMBS

COMBT=../data/total_combined
if [ -d $COMBT ]; then rm -r $COMBT; fi
mkdir -p $COMBT

METAPE=../meta/fastq_pe.txt
METASE=../meta/fastq_se.txt



# combine fastq files for total RNA libraries, renaming using namefunPE function
namefunPE () {
	OUT=${3}_${4}_${5}_${7}.fastq.gz
	echo $OUT
}

for i in {1..96}
do 
	F=$(sort -k 6,6 $METAPE | sed -n ${i}p | cut -f 1)

	OUT=$(namefunPE $(sort -k 6,6 $METAPE | sed -n ${i}p))

	echo $OUT $F
	
	cat $TOTAL/$F >>$COMBT/$OUT

done

# combine fastq files for small RNA libraries, renaming using namefunSE function
namefunSE () {
	OUT=${2}_${3}_${4}.fastq.gz
	echo $OUT
}

for i in {1..116}
do 
	F=$(sort -k 5,5 $METASE | sed -n ${i}p | cut -f 1)

	OUT=$(namefunSE $(sort -k 5,5 $METASE | sed -n ${i}p))

	echo $OUT $F
	
	cat $SMALL/$F >>$COMBS/$OUT

done

# some data files were uncompressed or already concatenated. 
# these are not in the metadata file and are dealt with separately here 
# deal with uncompressed files
cat $SMALL/*Hsa*fastq | gzip >$COMBS/Hsa_CONTROL_S1.fastq.gz
echo Hsa control done

cat $SMALL/*RE-8NoTxt* | gzip >$COMBS/RE_CONTROL_S7.fastq.gz
echo RE control done

gzip -c $SMALL/SNU_TruSeq13_R1.fastq >$COMBS/SNU_CONTROL_X.fastq.gz
echo SNU control done

