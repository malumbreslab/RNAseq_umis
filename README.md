# RNAseq_umis

## Description

RNA sequencing analysis pipeline with UMIs

Generate counts from fastq files from Ilumina sequencing with the use of UMIs.

## Table of contents

- [Workflow](#workflow)
- Contents of the repository

## Workflow

![This is an image](/images/workflow.png)

## Contents of the repository

- The bash script `rnaseq_pipeline_umis.sh` that can be used for obtaining read counts from fastq files.
- `Resource` folder has two files, `polyA.fa.gz` (extract polyA in trimmering) and `truseq_rna.fa.gz` (extract adapters in trimmering).
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
  
### Step 4: Create folders for result files

```
mkdir counts genome fastqc
```
### Step 5: Download genome reference and annotations files (hg38)

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
  
### Step 6: Index genome reference (hg38)

```
cd ..
STAR --runThreadN 16 --runMode genomeGenerate --genomeDir genome --genomeFastaFiles genome/genome_hg38.fa --sjdbGTFfile genome/annotations_hg18.gtf --sjdbOverhang 75
```

### Step 6a: Run pipeline in local

You must run this in your terminal shell and in sample must type the sample name, such as S1.

```
bash rnaseq_pipeline_umis.sh <sample>
```
  
### Step 6b: Run pipeline in cluster of CNIO

```
sbatch --mem=64G -t1440 -c 16 -o log.txt -e error.txt --wrap "bash rnaseq_pipeline_umis.sh <sample>"
```
  
## RECOMENDATIONS

The minimal requeriments are:

- Memory (mem): `64Gb`
- Number of cores (c): `16`
- Time (t): `1440 min`

For fastq with 10 millions of reads the time is about 40 min (without index genome reference, which could be around 2 hours)
  
When you index genome reference, you don't have to do it again if you want align some samples with the same reference genome.




