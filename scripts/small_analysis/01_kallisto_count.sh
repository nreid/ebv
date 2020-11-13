#!/bin/bash
#SBATCH --job-name=counts
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
## quantify expression with kallisto	##	
##########################################

module load kallisto/0.44.0

# input, output directories
INDEX=../../genome/small_transcripts.kallisto
INDIR=../../results/05a_small_trimmed
OUTDIR=../../results/small_analysis/kallisto_counts
mkdir -p $OUTDIR

# a bash array containing the sample IDs
LIST=($(ls $INDIR/*trim.fastq.gz | sed 's/.*trimmed\///' | sed 's/.trim.*//'))

# get one sample ID using the slurm array task ID
SAM=${LIST[$SLURM_ARRAY_TASK_ID]}

kallisto quant \
	-i $INDEX \
	-o $OUTDIR/${SAM} \
	-t 8 \
	--single \
	$INDIR/${SAM}.trim.fastq.gz

date 

