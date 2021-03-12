#!/bin/bash
#SBATCH --job-name=shortstack
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 20
#SBATCH --mem=50G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date


# load software
module load samtools/1.10
module load bowtie/1.1.2 

# shortstack does not have a systemwide install yet. using my home directory

# add RNAfold to path
	# this DOES have a system-wide install. switch to module. 
PATH=${PATH}:/home/CAM/nreid/bin/ViennaRNA-2.4.17/src/bin/

SS=~/bin/ShortStack/ShortStack

# input/output variables

INDIR=../../results/05a_small_trimmed/
FASTQ=$(find $INDIR -name "*fastq*" | sort | tr "\n" " ")
GENOME=../../genome/total_genome.fa

OUTDIR=../../results/small_analysis/shortstack_results

# run shortstack

$SS --bowtie_cores 20 --sort_mem 50G --readfile $FASTQ --genomefile $GENOME --outdir $OUTDIR

