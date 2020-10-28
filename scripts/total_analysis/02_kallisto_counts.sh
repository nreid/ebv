#!/bin/bash
#SBATCH --job-name=counts
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=30G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --array=[0-11]
##SBATCH --mail-type=ALL
##SBATCH --mail-user=
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err

hostname
date

##########################################
## kallisto quantification algorithm	##	
##########################################

module load kallisto/0.44.0

# input, output directories
INDEX=../../genome/total_transcripts.kallisto
INDIR=../../results/05b_total_trimmed
OUTDIR=../../results/total_analysis/counts
mkdir -p $OUTDIR

# a bash array containing the sample IDs
LIST=($(ls $INDIR/*trim_R1.fastq.gz | sed 's/.*trimmed\///' | sed 's/_trim.*//'))

# get one sample ID using the slurm array task ID
SAM=${LIST[$SLURM_ARRAY_TASK_ID]}

kallisto quant \
	-i $INDEX \
	-o $OUTDIR/${SAM} \
	-t 8 \
	$INDIR/${SAM}_trim_R1.fastq.gz $INDIR/${SAM}_trim_R2.fastq.gz

date 

