#!/bin/bash
#SBATCH --job-name=quantifier.pl
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

# identify miRNAs with mirdeep2.pl

# load software

module load mirdeep2/2.0.0.8

# input/output variables

INDIR=../../results/small_analysis/mirdeep2_all
OUTDIR=../../results/small_analysis/mirdeep2_all/quant
mkdir -p $OUTDIR

# run quantifier.pl

# breaking my usual rule and cd'ing to the output directory
	# b/c mirdeep2 does not offer output directory specification???

cd $OUTDIR

quantifier.pl \
-p ../ref_precursors.fa \
-m ../ref_miRNA.fa \
-r ../reads.fa

