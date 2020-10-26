#!/bin/bash
#SBATCH --job-name=trimmo_total
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --array=[0-11]%12


hostname
date

module load Trimmomatic/0.39

INDIR=../data/total_combined

OUTDIR=../results/05b_total_trimmed
mkdir -p $OUTDIR

ADAPT=../meta/TruSeq_total.fa

SAMARRAY=($(ls $INDIR | sed 's/_R..fastq.gz//' | sort | uniq))

SAM=${SAMARRAY[$SLURM_ARRAY_TASK_ID]}
java -jar $Trimmomatic PE -threads 4 \
        $INDIR/${SAM}_R1.fastq.gz \
        $INDIR/${SAM}_R2.fastq.gz \
        $OUTDIR/${SAM}_trim_R1.fastq.gz $OUTDIR/${SAM}_singles_R1.fastq.gz \
        $OUTDIR/${SAM}_trim_R2.fastq.gz $OUTDIR/${SAM}_singles_R2.fastq.gz \
        ILLUMINACLIP:$ADAPT:2:30:10 \
        SLIDINGWINDOW:4:25 \
        MINLEN:45 \
        CROP:75

date