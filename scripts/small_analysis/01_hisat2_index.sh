#!/bin/bash
#SBATCH --job-name=index
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=20G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
##SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

##########################################
## create index file					##	
##########################################

INDIR=../../genome/

module load hisat2/2.2.0

hisat2-build -p 16 $INDIR/total_transcripts.fa $INDIR/total_transcripts.hisat2

date 