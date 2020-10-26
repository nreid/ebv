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
