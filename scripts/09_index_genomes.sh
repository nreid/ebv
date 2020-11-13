#!/bin/bash
#SBATCH --job-name=index_genomes
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=80G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --array=[0-11]
##SBATCH --mail-type=ALL
##SBATCH --mail-user=
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err

# this array script indexes each reference genome using a few different aligners


INDIR=../../genome/


# total transcriptome
if [[ "$SLURM_ARRAY_TASK_ID" -eq 0 ]]; then
module load kallisto/0.44.0
kallisto index -i $INDIR/total_transcripts.kallisto $INDIR/total_transcripts.fa
module unload kallisto/0.44.0

module load hisat2/2.2.0
hisat2-build -p 16 $INDIR/total_transcripts.fa $INDIR/total_transcripts.hisat2
module unload hisat2/2.2.0

module load samtools/1.9
samtools index $INDIR/total_transcripts.fa
module unload samtools/1.9

module load bwa/0.7.15
bwa index -p $INDIR/total_transcripts.bwa $INDIR/total_transcripts.fa
module unload bwa/0.7.15

fi

# small transcriptome
if [[ "$SLURM_ARRAY_TASK_ID" -eq 1 ]]; then
module load kallisto/0.44.0
kallisto index -i $INDIR/small_transcripts.kallisto $INDIR/small_transcripts.fa
module unload kallisto/0.44.0

module load hisat2/2.2.0
hisat2-build -p 16 $INDIR/small_transcripts.fa $INDIR/small_transcripts.hisat2
module unload hisat2/2.2.0

module load samtools/1.9
samtools index $INDIR/small_transcripts.fa
module unload samtools/1.9

module load bwa/0.7.15
bwa index -p $INDIR/small_transcripts.bwa $INDIR/small_transcripts.fa
module unload bwa/0.7.15

fi

# total genome
if [[ "$SLURM_ARRAY_TASK_ID" -eq 2 ]]; then

module load hisat2/2.2.0
hisat2-build -p 16 $INDIR/total_genome.fa $INDIR/total_genome.hisat2
module unload hisat2/2.2.0

module load samtools/1.9
samtools index $INDIR/total_genome.fa
module unload samtools/1.9

module load bwa/0.7.15
bwa index -p $INDIR/total_genome.bwa $INDIR/total_genome.fa
module unload bwa/0.7.15

fi

# small genome
if [[ "$SLURM_ARRAY_TASK_ID" -eq 3 ]]; then

module load hisat2/2.2.0
hisat2-build -p 16 $INDIR/small_genome.fa $INDIR/small_genome.hisat2
module unload hisat2/2.2.0

module load samtools/1.9
samtools index $INDIR/small_genome.fa
module unload samtools/1.9

module load bwa/0.7.15
bwa index -p $INDIR/small_genome.bwa $INDIR/small_genome.fa
module unload bwa/0.7.15


fi