#!/bin/bash
#SBATCH --job-name=multiqc
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=10G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

echo `hostname`

# input dirs
INS=../results/06a_fastqc_small_trimmed
INT=../results/06b_fastqc_total_trimmed

# output dirs
OUTS=../results/07a_multiqc_small_trimmed
mkdir -p $OUTS

OUTT=../results/07b_multiqc_total_trimmed
mkdir -p $OUTT

module load MultiQC/1.9

# small RNA
multiqc --outdir $OUTS $INS

# total RNA
multiqc --outdir $OUTT $INT

