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
mkdir -p $COMBS

COMBT=../data/total_combined
mkdir -p $COMBT




