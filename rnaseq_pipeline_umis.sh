#!/bin/sh

## are necesary umi_tools, bbtools, fasqc, STAR, samtools, HTseq, rseqc
## problems between samtools and umi_tools packages
## are necesary fastqc, genome, counts folders
## requeriments: 64GB 16 cores (-c) 1440 time

start=$(date +%s)

for i in $@
do
    echo -e "\033[0;35m---- Analyzing $i sample ----"
    mkdir $i
    cat data/*${i}_*.fastq.gz > data/$i.fastq.gz
    echo "Total reads: $(($(zcat data/$i.fastq.gz|wc -l)/4|bc))"

    echo "Quality control with fastqc ..."
    fastqc -o fastqc/ -t 4 data/$i.fastq.gz
    echo "Quality control done"

    echo -e "\033[0;37mExtracting UMIs with umi_tools ..."
    #umi_tools extract --stdin=data/$i.fastq.gz --bc-pattern=NNNNNNNNNN --log=$i/$i.processed.log --stdout $i/$i.processed.fastq.gz
    umi_tools extract --random-seed 1234 --extract-method=regex --bc-pattern='(?P<umi_1>.{6})(?P<discard_1>TATA).*' -I data/$i.fastq.gz -S $i/$i.processed.fastq.gz --log=$i/$i.processed.log
    echo "Extracting UMIs done"
    echo "Reads after UMIs extraction: $(($(zcat $i/$i.processed.fastq.gz|wc -l)/4|bc))"

    echo "Trimmering with bbduk ..."
    bbduk.sh in=$i/$i.processed.fastq.gz out=$i/$i.clean.fastq.gz ref=resources/adapters.fa.gz k=13 ktrim=r useshortkmers=t mink=5 qtrim=r trimq=10 minlength=20
    echo "Trimmering done"
    echo "Reads after trimmering: $(($(zcat $i/$i.clean.fastq.gz|wc -l)/4|bc))"

    echo "Quality control with fastqc ..."
    fastqc -o fastqc/ -t 4 $i/$i.clean.fastq.gz
    echo "Quality control done"

    ## Genome reference Index
    ## only if the index is not done
    ## 64GB 16 cores (-c) 1440 time
    ## genome and gtf file download with wget from ncbi
    ## wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/latest/hg38.fa.gz
    ## wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz
    ## g38_RefSeq.bed must be downloaded
    ## wget https://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/hg38_RefSeq.bed.gz
    ## STAR --runThreadN 16 --runMode genomeGenerate --genomeDir genome --genomeFastaFiles genome/genome_hg38.fa --sjdbGTFfile genome/annotations_hg18.gtf --sjdbOverhang 75

    gzip -d $i/$i.clean.fastq.gz

    echo "Alignment with STAR ..."
    STAR --runThreadN 16 --genomeDir /storage/scratch01/shared/indexes/h_sapiens/gencode/v38/star/2.7.9a --readFilesIn $i/$i.clean.fastq --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.6 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMattributes NH HI NM MD --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $i/${i}.
    echo "Alignment done"
    echo "Mapped reads: $(samtools view -F 0x40 $i/$i.Aligned.sortedByCoord.out.bam | cut -f1 | sort | uniq | wc -l)"

    echo "Indexing with samtools ..."
    samtools index $i/$i.Aligned.sortedByCoord.out.bam
    echo "Indexing done"

    echo "Deduplication with umi_tools ..."
    umi_tools dedup -I $i/$i.Aligned.sortedByCoord.out.bam --output-stats=$i/deduplicated -S $i/$i.deduplicated.bam
    echo "Deduplication done"
    echo "Reads after deduplication: $(samtools view -F 0x40 $i/$i.deduplicated.bam | cut -f1 | sort | uniq | wc -l)"

    echo "Indexing with samtools ..."
    samtools index $i/$i.deduplicated.bam
    echo "Indexing done"

    echo "Counts with htseq ..."
    htseq-count -m union -s yes -f bam -i gene_name -r pos $i/$i.deduplicated.bam /storage/scratch01/shared/genomes/h_sapiens/gencode/v38/annotations.gtf.gz > counts/$i.read_counts.txt
    echo "Counts done"

    echo "Counts distribution with rseqc ..."
    read_distribution.py -i $i/$i.deduplicated.bam -r resources/hg38_RefSeq.bed > counts/$i.read_distribution.txt
    echo "Counts distribution done"

    echo "Sample $i done"
done

end=$(date +%s)
runtime=$((end-start))
echo "Execution time: $runtime seconds"
echo "Execution time: $((($runtime)/60)) minutes"

echo -e "\033[0;33mAll samples analyzed!"
