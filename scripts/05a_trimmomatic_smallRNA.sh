#!/bin/bash
#SBATCH --job-name=trimmo_small
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --array=[0-28]%28


hostname
date

module load Trimmomatic/0.39

INDIR=../data/small_combined

OUTDIR=../results/05a_small_trimmed
mkdir -p $OUTDIR

ADAPT=../meta/TruSeq_SmallRNA.fa

SAMARRAY=($(ls $INDIR | sed 's/.fastq.gz//'))

SAM=${SAMARRAY[$SLURM_ARRAY_TASK_ID]}

java -jar $Trimmomatic SE \
        -threads 4 \
        $INDIR/${SAM}.fastq.gz \
        $OUTDIR/${SAM}.trim.fastq.gz \
        ILLUMINACLIP:${ADAPT}:2:30:6 \
        SLIDINGWINDOW:4:20 \
        MINLEN:15 \
        CROP:75

date