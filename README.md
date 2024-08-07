# RNAseq_umis

## Description

Pipeline to analyze 3' mRNA sequencing from FFPE samples and using UMIs

## Table of contents

- [Workflow](#workflow)
- [Contents of the repository](#contents-of-the-repository)
- [Pipeline](#pipeline)
- [Recomendations](#recomendations)

## Workflow

![This is an image](/images/workflow.png)

## Contents of the repository

- The bash script `rnaseq_pipeline_umis.sh` that can be used for obtaining read counts from fastq files.
- `Resource` folder has `adapters.fa.gz` with adapters sequences and polyA.
- `config_env.yml` which is a file for create the working enviroment.

The output of this bash script includes:

- `<sample>.read_counts.txt` which are the read counts (counts folder).
- `<sample>.read_distribution.txt` which are the read counts (counts folder).
- `<sample>.deduplicated.bam` which are deduplicated bam (<sample> folder).

There are other outputs less important such as fastqc, bai, trimmering, umi data, ...

## Pipeline
  
### Step 1: Clone the repo in your home

If you have not wget package installed
  
```
conda install -c anaconda git
```
  
Then clone the entire repository in your local space

```
git clone https://github.com/malumbreslab/RNAseq_umis.git
cd RNAseq_umis
```
  
### Step 2: Update conda and create the environment

```
conda update --all
conda env create -f config_env.yml
```
  
### Step 3: Activate the environment

```
conda activate rseq
```
### Step 4: Create folders for data

In this folder you must introduce compressed fastq files (Example: S1.fastq.gz)
  
```
mkdir data
```
  
### Step 5: Create folders for result files

```
mkdir counts genome fastqc
```
### Step 6: Download genome reference and annotations files (hg38)

Steps 6 and 7 are not neccesary if have a created index genome (or use the index genome from shared folder in the cluster)

If you have not wget package installed
  
```
conda install -c anaconda wget
```
Then download files with url and uncompress
  
```
cd genome
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/latest/hg38.fa.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz
wget https://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/hg38_RefSeq.bed.gz  
gzip -d hg38.fa.gz
gzip -d hg38.ncbiRefSeq.gtf.gz
gzip -d hg38_RefSeq.bed.gz 
```
  
### Step 7: Index genome reference (hg38)

```
cd ..
STAR --runThreadN 16 --runMode genomeGenerate --genomeDir genome --genomeFastaFiles genome/genome_hg38.fa --sjdbGTFfile genome/annotations_hg18.gtf --sjdbOverhang 75
```

### Step 8a: Run pipeline in local

You must run this in your terminal shell and in sample must type after the script the diferent samples names that you want to analyze separated by spaces.

```
bash rnaseq_pipeline_umis.sh <sample1> <sample2> <sample3>
```
  
### Step 8b: Run pipeline in cluster of CNIO

Outputs:
  
- `log.txt` is the output file
- `error.txt` is the error file
  
Parameters:
  
- `--mem` is memory
- `-t` is time
- `-J` is job name
- `-c` is number of cores
- `-o` is name of output file
- `-e` is name of error file
- `--wrap` is the command which you want to run on the cluster
  
```
sbatch --mem=64G -t1440 -c 16 -J name -o log.txt -e error.txt --wrap "bash rnaseq_pipeline_umis.sh <sample1> <sample2> <sample3>"
```

## Recomendations

The minimal requeriments are:

- Memory (mem): `64Gb`
- Number of cores (c): `16`
- Time (t): `1440 min`

For fastq with 10 millions of reads the time is about 40 min (without index genome reference, which could be around 2 hours)
  
When you index genome reference, you don't have to do it again if you want align some samples with the same reference genome.

## Important functions

### Quality control

It is used fastqc

### UMIs

umi_tools (extract/dedup) is used in order to extract and deduplicate reads

### Trimmering

bbduk.sh is used in order to remove adapters and polyA tail from reads

### Aligment

It is recommended use STAR for RNA aligment

### Counting

Htseq-count is used -m union (more reads), -s yes (stranded reads), -i gene_name (hugo symbol)








