#!/bin/bash 
#SBATCH --job-name=symlink_data
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


# symlink the raw data to this github repo

DATADIR=../data
SMALL=$DATADIR/small_fastq
TOTAL=$DATADIR/total_fastq

mkdir -p $SMALL
mkdir -p $TOTAL

# total rna libraries
RAWDIR1=/labs/Oneill/jules/paired_end_fastas/
ln -s $RAWDIR1/*fastq.gz $TOTAL

# single end data
RAWDIR2=/labs/Oneill/jules/single_end_fastas/
ln -s $RAWDIR2/*fastq.gz $SMALL

# several single end files are uncompressed for some reason. 
ln -s $RAWDIR2/*fastq $SMALL

######## NOTE:
######## THERE IS A TYPO IN ONE FILE NAME:
	####### RE-8-TruSeq_DroskaKD_S1_L003_R1_001.fastq.gz 
	####### SHOULD BE:
	####### RE-8-TruSeq_DroshaKD_S1_L003_R1_001.fastq.gz
	####### to correct:  Droska -> Drosha

	####### this error cannot have been made by the sequencer, right?


	####### this is dealt with downstream by combining fastq files per the metadata tables. 

