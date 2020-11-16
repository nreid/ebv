#!/bin/bash
#SBATCH --job-name=bwa_map
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=30G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --array=[0-28]
##SBATCH --mail-type=ALL
##SBATCH --mail-user=
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err

hostname
date

##########################################
## quantify expression with bwa	##	
##########################################

module load bwa/0.7.17
module load samtools/1.9

# input, output directories
INDEX=../../genome/small_transcripts.bwa
INDIR=../../results/05a_small_trimmed
OUTDIR=../../results/small_analysis/bwa_trans_map
mkdir -p $OUTDIR

# a bash array containing the sample IDs
LIST=($(ls $INDIR/*trim.fastq.gz | sed 's/.*trimmed\///' | sed 's/.trim.*//'))

# get one sample ID using the slurm array task ID
SAM=${LIST[$SLURM_ARRAY_TASK_ID]}

bwa aln $INDEX $INDIR/${SAM}.trim.fastq.gz | \
samtools view -S -h -u - | \
samtools sort -T $SAM - >$OUTDIR/${SAM}.bam

samtools index $OUTDIR/${SAM}.bam

date 

