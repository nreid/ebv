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
INS=../results/03a_fastqc_small
INT=../results/03b_fastqc_total

# output dirs
OUTS=../results/04a_multiqc_small
mkdir -p $OUTS

OUTT=../results/04b_multiqc_total
mkdir -p $OUTT

module load MultiQC/1.9

# small RNA
multiqc --outdir $OUTS $INS

# total RNA
multiqc --outdir $OUTT $INT

