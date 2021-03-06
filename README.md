# ebv

## Background

original data located here
`/labs/Oneill/jules`

new project directory here:
`/labs/CBC/core_projects/2020/ebv/`


There are two datasets, a set of total RNA libraries and a set of small RNA libraries. They will be analyzed in parallel, at least at the beginning. 

adapters for trimming obtained [on illumina's website](https://support.illumina.com/bulletins/2016/12/what-sequences-do-i-use-for-adapter-trimming.html)
total rna libraries use the TruSeq adapters
small rna libraries use TruSeq small RNA adapter

To reproduce this analysis on the UConn Xanadu cluster, ensure the raw data is still there, then run each script in sequence. 

## Analysis steps

### Step 1: Obtain data

The data are located in this directory on xanadu: `/labs/Oneill/jules`. [`01_symlinkrawdata.sh`](/scripts/01_symlinkrawdata.sh) symlinks the data to the `/data` directory. 

### Step 2: Combine data across lanes

To make downstream analysis easier, I combine sequencing lanes within samples. This is accomplished using the script [`02_combine_lanes.sh`](/scripts/02_combine_lanes.sh) and metadata tables in the directory [meta](/meta). Inconsistent file name formatting is fixed here, along with one sample name typo. 

### Step 3: Run `FastQC`

To assess sequence quality, I run `fastqc` on both sets of libraries with the scripts [`03a_fastqc_small.sh`](/scripts/03a_fastqc_small.sh) and [`03b_fastqc_total.sh`](/scripts/03b_fastqc_total.sh). 

### Step 4: Run `MultiQC`

I aggregate `fastqc` reports with `multiqc` using [`04_multiqc.sh`](/scripts/04_multiqc.sh)

### Step 5: Run `Trimmomatic`

To remove adapter contamination and low quality sequence, I run `Trimmomatic` on both sets of libraries, with adapter fasta files located in [`meta`](/meta) and scripts [`05a_trimmomatic_smallRNA.sh`](/scripts/05a_trimmomatic_smallRNA.sh) and [`05b_trimmomatic_totalRNA.sh`](/scripts/05b_trimmomatic_totalRNA.sh). 

### Step 6: Run `FastQC` again

To assess trimmed sequence quality, I run `fastqc` on both sets of libraries with the scripts [`06a_fastqc_small_trimmed.sh`](/scripts/06a_fastqc_small_trimmed.sh) and [`06b_fastqc_total_trimmed.sh`](/scripts/06b_fastqc_total_trimmed.sh). 

### Step 7: Run `MultiQC` again

I aggregate trimmed `fastqc` reports with `multiqc` using [`07_multiqc_trimmed.sh`](/scripts/07_multiqc_trimmed.sh)

### Step 8: Get and format genomes and annotations

We need to aggregate and modify reference files to map against. The total RNA libraries have ERCC spike-ins added. The small RNA libraries have ExiSeq small RNA spike-ins added. All samples (except controls) are either infected by, or transfected with, Epstein-Barr virus in full or in part. 

I use the script [`08_get_genomes.sh`](/scripts/08_get_genomes.sh). 

In short, I pull a bunch of files from NCBI, ENSEMBL, and the manufacturer of the spike-ins (with the exception of the ExiSEQ sequences, which I couldn't find, so I'm using the file from `/labs/Oneill/jules`). 

I generate modified reference transcriptomes by appending the ENSEMBL transcripts for coding genes and ncRNA with the NCBI EBV transcripts and either the small RNA spike-ins for the small RNA analysis (`small_transcripts.fa`) or the ERCC spike-ins for the total RNA analysis (`total_transcripts.fa`). 

This script will be updated to create a modified reference genome/annotation pair as well to do genome-based quantification. 

______________


The analysis pathways now diverge sufficiently that they will be documented separately. 

## Total RNA libraries

First I will quantify gene expression on the synthetic reference transcriptome created above (`total_transcripts.fa`) using `kallisto` and analyzing the resulting count data in `DESeq2`. 

### Step 1: Index the reference file

See script [`01_kallisto_index.sh`](/scripts/total_analysis/01_kallisto_index.sh)

### Step 2: Quantify expression

See script [`02_kallisto_counts.sh`](scripts/total_analysis/02_kallisto_counts.sh)

### Step 3: Statistical analysis in R

See script [`analysis.R`](scripts/total_analysis/analysis.R)

## Small RNA libraries

I have been exploring a few different approaches here:

`ShortStack`: See script [`01_shortstack.sh`](https://github.com/nreid/ebv/blob/main/scripts/small_analysis/01_shortstack.sh). This script runs [`ShortStack`](https://github.com/MikeAxtell/ShortStack) on all the small RNA samples. It appears to have run successfully. 

`miRDeep2`: [`mirDeep2`](https://github.com/rajewsky-lab/mirdeep2) is a a series of scripts. 
 - [`02a_mirdeep2_bowtie_build.sh`](https://github.com/nreid/ebv/blob/main/scripts/small_analysis/02a_mirdeep2_bowtie_build.sh) This step indexes the reference genome (+EBV). 
 - [`02b_mirdeep2_mapper.pl.sh`](https://github.com/nreid/ebv/blob/main/scripts/small_analysis/02b_mirdeep2_mapper.pl.sh) This step preprocesses the reads and maps them to the reference genome. it produces a file of unique sequences (`reads.fa`) and a mapping of those reads to the genome (`reads_vs_genome.arf`). 
 - [`02c_mirdeep2_mirdeep2.pl.sh`](https://github.com/nreid/ebv/blob/main/scripts/small_analysis/02c_mirdeep2_mirdeep2.pl.sh) This is the core script, it takes in a file of known miRNAs (and their precursors, and miRNAs from related species), decides whether they are present in the current dataset and searches for new candidate miRNAs. 
 - [`02d_mirdeep2_quantifier.pl.sh`](https://github.com/nreid/ebv/blob/main/scripts/small_analysis/02d_mirdeep2_quantifier.pl.sh) I used this script to quantify miRNAs of interest. `mirdeep2.pl` quantifies miRNAS, but for some reason it only output the human-specific counts, and did not include the EBV counts. 

I have also aligned all the small RNA libraries to the reference genome (including EBV) with `bwa`. That is documented here: [`xx_bwa_genome_map.sh`](https://github.com/nreid/ebv/blob/main/scripts/small_analysis/xx_bwa_genome_map.sh). 

PCA and heatmap plots and a description of some results from both of these datasets can be found [here](). 

