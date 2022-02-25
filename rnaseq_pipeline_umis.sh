#!/bin/sh

## are necesary umi_tools, bbtools, fasqc, STAR, samtools, HTseq, rseqc
## problems between samtools and umi_tools packages
## are necesary fastqc, genome, counts folders (mkdir fasqc, genome, counts)
## requeriments: 64GB 16 cores (-c) 1440 time

echo 'All samples in one ...'
mkdir $1
cat data/*$1_*.fastq.gz > data/$1.fastq.gz
echo 'All samples in one done'

## only if you want uncompress the fasta file
#gzip -d data/$1.fastq.gz

echo 'Extracting UMIs with umi_tools ...'
umi_tools extract --stdin=data/$1.fastq.gz --bc-pattern=NNNNNNNNNN --log=$1/$1.processed.log --stdout $1/$1.processed.fastq.gz
echo 'Extracting UMIs done'

echo 'Trimmering with bbduk ...'
bbduk.sh in=$1/$1.processed.fastq.gz out=$1/$1.clean.fastq.gz ref=resources/polyA.fa.gz,resources/truseq_rna.fa.gz k=13 ktrim=r useshortkmers=t mink=5 qtrim=r trimq=10 minlength=20
echo 'Trimmering done'

echo 'Quality control with fasqc...'
fastqc -o $1/ -t 4 $1/$1.clean.fastq.gz
echo 'Quality control done'

## only if the index is not done
## 64GB 16 cores (-c) 1440 time
## genome and gtf file download with wget from ncbi
## wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/latest/hg38.fa.gz
## wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz
## STAR --runThreadN 16 --runMode genomeGenerate --genomeDir genome --genomeFastaFiles genome/genome_hg38.fa --sjdbGTFfile genome/annotations_hg18.gtf --sjdbOverhang 75

echo 'Uncompressing file with gzip...'
gzip -d $1/$1.clean.fastq
echo 'Uncompressing file done'

echo 'Alignment with STAR...'
STAR --runThreadN 16 --genomeDir genome --readFilesIn $1/$1.clean.fastq --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.6 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMattributes NH HI NM MD --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $1/$1
echo 'Alignment done'

echo 'Indexing with samtools...'
samtools index $1/$1Aligned.sortedByCoord.out.bam
echo 'Indexing done'

echo 'Deduplication with umi_tools...'
umi_tools dedup -I $1/$1Aligned.sortedByCoord.out.bam --output-stats=$1/deduplicated -S $1/$1.deduplicated.bam
echo 'Deduplication done'

echo 'Counts with htseq...'
htseq-count -m intersection-nonempty -s yes -f bam -r pos $1/$1.deduplicated.bam genome/annotations_hg18.gtf > counts/$1.read_counts.txt
echo 'Counts done'

## g38_RefSeq.bed must be downloaded
## wget https://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/hg38_RefSeq.bed.gz
echo 'Counts distribution with rseqc...'
read_distribution.py -i $1/$1.deduplicated.bam -r genome/hg38_RefSeq.bed > counts/$1.read_distribution.txt
echo 'Counts distribution done'

## Total counts
awk '{s+=$2}END{print s}' counts/$1.read_counts.txt