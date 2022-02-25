# RNAseq_umis

RNA sequencing analysis pipeline

## Prerequisites

First of all you must have conda installed and you have to create a new enviroment (using the config_env.yml file).

`conda env create -f config_env.yml`

`conda activate rseq`

## Using pipeline

Create directory of working and enter.

`mkdir rnaseq`

`cd rnaseq`

Create folder for save result files

`mkdir fasqc genome counts`

Download genome reference and annotations files (this url is for hg38 genome reference)

`cd genome`

`wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/latest/hg38.fa.gz`

`wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz`

`wget https://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/hg38_RefSeq.bed.gz`

Uncompress files

`gzip -d hg38.fa.gz`

`gzip -d hg38.ncbiRefSeq.gtf.gz`

`gzip -d hg38_RefSeq.bed.gz`

Index genome reference

`cd ..`

`STAR --runThreadN 16 --runMode genomeGenerate --genomeDir genome --genomeFastaFiles genome/genome_hg38.fa --sjdbGTFfile genome/annotations_hg18.gtf --sjdbOverhang 75`







