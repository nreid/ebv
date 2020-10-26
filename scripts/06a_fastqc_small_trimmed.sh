#!/bin/bash
#SBATCH --job-name=fastqc_small
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --array=[0-28]%29

hostname
date

echo "host name : " `hostname`
echo This is array task number $SLURM_ARRAY_TASK_ID

# load software
module load fastqc/0.11.7

#input/output directories, supplementary files

FASTQS=($(find ../results/05a_small_trimmed -name "*fastq.gz"))

OUTDIR=../results/06a_fastqc_small_trimmed
mkdir -p $OUTDIR

INFILE=${FASTQS[$SLURM_ARRAY_TASK_ID]}

# run fastqc. "*fq" tells it to run on all fastq files in directory "../rawdata/"
fastqc -t 2 -o $OUTDIR $INFILE

date