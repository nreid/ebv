#!/bin/bash
#SBATCH --job-name=bowtie_build_mirdeep2
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

# build a bowtie index for mirdeep2

module load mirdeep2/0.1.3

GENOME=../../genome/total_genome.fa

OUT=../../genome/total_genome_mirdeep2

bowtie-build $GENOME $OUT