#!/bin/bash
#SBATCH --job-name=get_genomes
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=10G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

# GOALS: 
	# 1. a human genome fasta with ERCC and ExiSEQ spike-ins
	# 2. corresponding annotation files in GFF3 and GTF format
	# 3. a ncRNA transcript file with ExiSEQ spike-ins
	# 4. a CDS transcript file with ERCC spike-ins


# We need a few different reference files, and they need to be modified 
# so that we can map to spike-ins and the epstein-barr virus. 

# We will use human genome resources from here:
	# human genome GRCh38.p13
	# ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/GRCh38_major_release_seqs_for_alignment_pipelines

	# we'll use the no_alts genome and download the corresponding indexes: 
		# HISAT2
		# BWA
		# samtools

	# we'll also download the GTF and GFF3 annotation files

	# THIS HUMAN GENOME ALREADY CONTAINS A COPY OF THE EPSTEIN BARR VIRUS, "chrEBV"
		# the annotation file does not contain EBV annotations though. 
		# we need to obtain the annotation and append it. 
		# this is the reference assembly for chrEBV:
			# https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/402/265/GCF_002402265.1_ASM240226v1/
			# we will validate it before appending the annotations. 

# as an alternate to the genome-based approach, we'll use ensembl transcripts as well. 
	# https://uswest.ensembl.org/info/data/ftp/index.html



GENDIR=../genome
mkdir -p $GENDIR

cd $GENDIR

NCBIHOM=ncbi_homo
NCBIEBV=ncbi_ebv
ENSEMBL=ensembl_homo
SPIKE=spike_ins
mkdir -p $NCBIHOM $NCBIEBV $ENSEMBL $SPIKE


############################
# Get human genome files
# ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/GRCh38_major_release_seqs_for_alignment_pipelines
############################

cd $NCBIHOM 

# get reference genome:
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

# get annotations in gff3 and gtf format:
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf.gz
gunzip GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gff.gz
gunzip GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gff.gz

# get samtools index
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai

# # get hisat2 index
# wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.hisat2_index.tar.gz
# tar -xvf GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.hisat2_index.tar.gz

# # get bwa index
# wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bwa_index.tar.gz
# tar -xvf GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bwa_index.tar.gz

############################
# Get EBV genes and annotation
# https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/402/265/GCF_002402265.1_ASM240226v1/
############################

cd ../$NCBIEBV

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/402/265/GCF_002402265.1_ASM240226v1/GCF_002402265.1_ASM240226v1_genomic.gff.gz
gunzip GCF_002402265.1_ASM240226v1_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/402/265/GCF_002402265.1_ASM240226v1/GCF_002402265.1_ASM240226v1_genomic.gtf.gz
gunzip GCF_002402265.1_ASM240226v1_genomic.gtf.gz

# RNAs, including small RNAs and large polycistrons
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/402/265/GCF_002402265.1_ASM240226v1/GCF_002402265.1_ASM240226v1_rna_from_genomic.fna.gz

############################
# Get ENSEMBL transcripts
# https://uswest.ensembl.org/info/data/ftp/index.html
############################

cd ../$ENSEMBL

# cdna
wget ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
gunzip Homo_sapiens.GRCh38.cdna.all.fa.gz
# ncRNA
wget ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz
gunzip Homo_sapiens.GRCh38.ncrna.fa.gz

############################
# Get spike in sequences and annotations
############################

cd ../$SPIKE

# these are found in /labs/Oneill/jules/jules_stuff
	# they don't have good provenance

# small RNA spike-ins
cp /labs/Oneill/jules/jules_stuff/ExiSEQ-NGS-QC-Spike-ins.fa .

# ercc spike-ins. 
	# https://www.thermofisher.com/order/catalog/product/4456740?us&en#/4456740?us&en
	# readme
	wget https://assets.thermofisher.com/TFS-Assets/LSG/manuals/cms_095048.txt
	# concentrations in the mixes
	wget https://assets.thermofisher.com/TFS-Assets/LSG/manuals/cms_095046.txt
	# sequences
	wget https://assets.thermofisher.com/TFS-Assets/LSG/manuals/cms_095047.txt

	awk '{OFS=""}{print ">",$1,"\n",$5}' cms_095047.txt >ercc_spike_ins.fa

### drop back to GENDIR

cd ..

############################
# Create modified human genome and annotation
############################

# Add ERCC spike-ins, ExiSEQ spike-ins to GRCh38.p13 (EBV is already there)

# add spike in lines to GTF and GFF3 annotations




############################
# Create modified human transcriptomes
############################

# Add ERCC spike-ins and ncRNA to Homo_sapiens.GRCh38.cdna.all.fa.gz for total RNA analysis

cat \
$ENSEMBL/Homo_sapiens.GRCh38.cdna.all.fa \
$ENSEMBL/Homo_sapiens.GRCh38.ncrna.fa \
$SPIKE/ercc_spike_ins.fa \
>total_transcripts.fa


# Add ExiSEQ spike-ins to Homo_sapiens.GRCh38.ncrna.fa.gz for small RNA analysis

cat \
$ENSEMBL/Homo_sapiens.GRCh38.cdna.all.fa \
$ENSEMBL/Homo_sapiens.GRCh38.ncrna.fa \
$SPIKE/ExiSEQ-NGS-QC-Spike-ins.fa \
>small_transcripts.fa
