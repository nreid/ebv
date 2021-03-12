#!/bin/bash
#SBATCH --job-name=mirdeep2_mapper.pl
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 32
#SBATCH --mem=50G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

# map reads from all samples to reference genome with mapper.pl

# load software

module load mirdeep2/0.1.3

MAPPER=/isg/shared/apps/mirdeep2/0.1.3/src/mapper.pl

# input/output variables

INDIR=../../results/05a_small_trimmed/
OUTDIR=../../results/small_analysis/mirdeep2_all
FQDIR=$OUTDIR/fastq
mkdir -p $FQDIR

BOWTIE_GENOME_INDEX=../../genome/total_genome_mirdeep2

# mirdeep2 appears to want only decompressed fq files

for file in $(find $INDIR -name "*fastq*" | sort); do echo gunzip -c -d $file \>$FQDIR/$(echo $file | basename $file .gz); done | \
parallel -j 32

# for multiple samples, need a config file formatted like:

# sequencing_data_sample1.fa  sd1
# sequencing_data_sample2.fa  sd2
# sequencing_data_sample3.fa  sd3

# this line makes this file for all 32 samples
paste \
<(find $FQDIR -name "*fastq*" | sort) \
<(echo {001..032} | tr " " "\n") \
>$OUTDIR/config.txt

# run mirdeep2 mapper.pl script

# -d, -e mean config file provided, fastq files
# -i, -j convert from rna to dna, remove seqs with non-IUPAC characters
# -l discard reads < 18
# -m collapse reads
# -p map to genome, provide index
perl $MAPPER $OUTDIR/config.txt -d -e -h -i -j -l 18 -m -p $BOWTIE_GENOME_INDEX -s $OUTDIR/reads.fa -t $OUTDIR/reads_vs_genome.arf 

