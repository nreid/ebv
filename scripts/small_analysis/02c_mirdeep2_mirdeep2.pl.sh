#!/bin/bash
#SBATCH --job-name=mirdeep2.pl
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

# identify miRNAs

# load software

module load mirdeep2/0.1.3

MIRDEEP2=/isg/shared/apps/mirdeep2/0.1.3/src/miRDeep2.pl

# input/output variables

INDIR=../../results/small_analysis/mirdeep2_all

GENOME=../../genome/total_genome.fa

INFA=$INDIR/reads.fa
INARF=$INDIR/reads_vs_genome.arf


# get the miRNA database from mirbase
# extract human (hsa) and ebv (ebv) mirnas and their precursors 
# also extract pan (ptr) mirnas

wget -P $INDIR ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz 
	# homo/ebv
	zcat $INDIR/mature.fa.gz | grep -A 1 ">hsa" | grep -v "\-\-" | sed 's/ .*// ' >$INDIR/ref_miRNA.fa
	zcat $INDIR/mature.fa.gz | grep -A 1 ">ebv" | grep -v "\-\-" | sed 's/ .*// ' >>$INDIR/ref_miRNA.fa
	# pan 
	zcat $INDIR/mature.fa.gz | grep -A 1 ">ptr" | grep -v "\-\-" | sed 's/ .*// ' >$INDIR/ptr_miRNA.fa


wget -P $INDIR ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.gz
	# homo/ebv precursors
	zcat $INDIR/hairpin.fa.gz | grep -A 1 ">hsa" | grep -v "\-\-" | sed 's/ .*// ' >$INDIR/ref_precursors.fa
	zcat $INDIR/hairpin.fa.gz | grep -A 1 ">ebv" | grep -v "\-\-" | sed 's/ .*// ' >>$INDIR/ref_precursors.fa

# mirdeep2 cannot handle a genome fa with whitespace in the headers
# so we have to make another copy of the genome

GENOME2=../../genome/total_genome_nowhitespace.fa

sed 's/ .*//' $GENOME >$GENOME2

# run mirdeep2.pl

perl $MIRDEEP2 \
$INFA \
$GENOME2 \
$INARF \
$INDIR/ref_miRNA.fa \
$INDIR/ptr_miRNA.fa \
$INDIR/ref_precursors.fa \
-P \
-v \
-t Human \
2>$INDIR/report.log 

