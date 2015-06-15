#!/bin/bash

if [ -z "$1" ]; then 
	echo usage: $0 single_directory_of_raw_reads genome_build
    exit
fi

WD=$1
Build=$2
echo Genome Build Set to: $Build
echo Processing files in : $WD

cd $WD || { echo ERROR: could not find $WD , exiting... ; exit 1; }
mkdir raw || { echo 'ERROR: could not make raw directory, exiting...' ; exit 1; }
gunzip -v *.gz || { echo 'ERROR: could not unzip .fastq.gz files, exiting...' ; exit 1; }
echo Concatenating R1...
cat *_R1_0??.fastq > ./R1.fastq || { echo 'ERROR: could not concatenate reads into one file, exiting...' ; exit 1; }

echo Trimming Off First 3 Bases...
/home/clf21/bin/fastx_trimmer -Q33 -f 4 -l 50 -i ./R1.fastq -o ./R1_3Trimmed.fastq || { echo 'ERROR: could not trim reads, exiting...' ; exit 1; }

echo Trimming out Illumina Truseq Adapter Sequences...
java -classpath /home/clf21/bin/Trimmomatic-0.22/trimmomatic-0.22.jar org.usadellab.trimmomatic.TrimmomaticSE -threads 6 -phred33 ./R1_3Trimmed.fastq R1_filtered.fastq ILLUMINACLIP:/home/clf21/bin/TruSeqAdapters.fa:2:40:15 MINLEN:46

echo Producing QC report on raw fastq file...
/home/clf21/FastQC/fastqc --noextract ./R1.fastq
echo Producing QC report after trimming and filtering out adapters...
/home/clf21/FastQC/fastqc --noextract ./R1_filtered.fastq

echo Cleaning Up...
gzip *R?_0??.fastq
mv *R?_0??.fastq.gz ./raw/
mv *.csv ./raw/
rm ./R?.fastq
rm R1_3Trimmed.fastq
mkdir tophat || { echo 'ERROR: could not make tophat directory for output, exiting...' ; exit 1; }
mkdir logs || { echo 'ERROR: could not make output directory for program logs, exiting...' ; exit 1; }

echo Beginning Alignment...
/home/clf21/samtools/bin/tophat --segment-length 23 --library-type fr-unstranded -p 4 --transcriptome-index=/home/clf21/RNA-seq/GTF_knownGenes/compiled_$Build -x 4 -n 2 -o ./tophat/ /home/clf21/bin/bowtie-0.12.7/indexes/$Build ./R1_filtered.fastq 2> ./logs/tophat_log.txt
# this is for NO novel transcripts: /home/clf21/samtools/bin/tophat --segment-length 23 --library-type fr-unstranded -p 4 --transcriptome-index=/home/clf21/RNA-seq/GTF_knownGenes/compiled_$Build -T -x 4 -n 2 -o ./tophat/ /home/clf21/bin/bowtie-0.12.7/indexes/$Build ./R1_filtered.fastq 2> ./logs/tophat_log.txt

echo Cleaning Up...
rm ./R1_filtered.fastq
mkdir cufflinks  || { echo 'ERROR: could not make cufflinks directory for output, exiting...' ; exit 1; }

echo Beginning Cufflinks FPKM Estimates...
echo Using Reference Transcriptome found at /home/clf21/RNA-seq/GTF_knownGenes/"$Build"_genes.gtf ...
/home/clf21/bin/cufflinks -p 4 -u -o ./cufflinks/ --GTF-guide /home/clf21/RNA-seq/GTF_knownGenes/"$Build"_genes.gtf ./tophat/accepted_hits.bam 2> ./logs/cufflinks_log.txt
# this is for NO novel transcripts: /home/clf21/bin/cufflinks -p 4 -u -o ./cufflinks/ -G /home/clf21/RNA-seq/GTF_knownGenes/"$Build"_genes.gtf ./tophat/accepted_hits.bam 2> ./logs/cufflinks_log.txt

echo Making bigWig coverage track for browser visualization...
/home/clf21/bin/genomeCoverageBed -split -bg -ibam ./tophat/accepted_hits.bam -g /home/clf21/bin/chrom.sizes."$Build".txt > ./tophat/accepted_hits.bedGraph
/home/clf21/bin/bedGraphToBigWig ./tophat/accepted_hits.bedGraph /home/clf21/bin/chrom.sizes."$Build".txt ./tophat/accepted_hits.bigWig
rm ./tophat/accepted_hits.bedGraph
echo Normalizing bigWig signal for comparable browser coverage tracks...
/home/clf21/bin/java_genomics_toolkit/toolRunner.sh wigmath.Scale -i ./tophat/accepted_hits.bigWig -o ./tophat/accepted_hits_norm.wig 2> ./logs/wigScale_log.txt
# By default, this will scale the bigWig coverage to a mean of 1
/home/clf21/bin/wigToBigWig ./tophat/accepted_hits_norm.wig /home/clf21/bin/chrom.sizes."$Build".txt ./tophat/accepted_hits_scaled.bigWig
rm ./tophat/accepted_hits_norm.wig

echo 'Done. Run appears to have completed successfully, but check for error messages from tophat and cufflinks.' 




