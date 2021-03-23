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

# identify miRNAs with mirdeep2.pl

# load software

module load mirdeep2/2.0.0.8

MIRDEEP2=/isg/shared/apps/mirdeep2/2.0.0.8/src/miRDeep2.pl

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
# it also can't handle IUPAC ambiguities. so we'll remove those as well. 

#sed '/^[^>]/ s/[^AGTC]/N/gi' < seq.fa

GENOME2=../../genome/total_genome_nowhitespace.fa

# first sed call removes white space
# second sed call edits non header lines and replaces non ACGT characters with N and is case-insensitive
sed 's/ .*//' $GENOME | sed '/^[^>]/ s/[^AGTC]/N/gi' >$GENOME2

# run mirdeep2.pl

# breaking my usual rule and cd'ing to the mirdeep directory
	# b/c mirdeep2.pl does not offer output directory specification???

cd $INDIR

perl $MIRDEEP2 \
reads.fa \
../../../genome/total_genome_nowhitespace.fa  \
reads_vs_genome.arf \
ref_miRNA.fa \
ptr_miRNA.fa \
ref_precursors.fa \
-P \
-t Human 


# the file "result*csv" contains 4 different tab separated results tables
	# here we're going to just extract the table of known miRNAs that were detected in the analysis
	# we ignore novel ones because none were found in EBV, which is the focus here

awk '{if ($1 ~ /mature/) x=1}
	{if ($1 ~ /^$/) x=0}
	{if (x == 1) print $0}' \
	result_*.csv \
>mature_miRNA_found.txt

