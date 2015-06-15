#!/bin/bash

if [ -z "$1" ]; then 
	echo usage: $0 '<single_directory_of_raw_reads> <genome_build> <trim_length>'
    exit
fi

WD=$1
Build=$2
TrimLength=$3
echo Genome Build Set to: $Build
echo Processing files in : $WD
echo Trimming off $TrimLength bases ...
TrimLengthPlusOne=$((TrimLength + 1))

cd $WD || { echo ERROR: could not find $WD , exiting... ; exit 1; }
mkdir raw || { echo 'ERROR: could not make raw directory, exiting...' ; exit 1; }
gunzip -v *.gz || { echo 'ERROR: could not unzip .fastq.gz files, exiting...' ; exit 1; }
echo Concatenating R1...
cat *_R1_0??.fastq > ./R1.fastq || { echo 'ERROR: could not concatenate Left reads into one file, exiting...' ; exit 1; }
echo Concatenating R2...
cat *_R2_0??.fastq > ./R2.fastq || { echo 'ERROR: could not concatenate Right reads into one file, exiting...' ; exit 1; }

echo Trimming Off First $TrimLength Bases...
fastx_trimmer -Q33 -f $TrimLengthPlusOne -l 50 -i ./R1.fastq -o ./R1_Trimmed.fastq || { echo 'ERROR: could not trim Left reads, exiting...' ; exit 1; }
fastx_trimmer -Q33 -f $TrimLengthPlusOne -l 50 -i ./R2.fastq -o ./R2_Trimmed.fastq || { echo 'ERROR: could not trim Right reads, exiting...' ; exit 1; }

echo Trimming out Illumina Truseq Adapter Sequences...
java -classpath /home/clf21/bin/Trimmomatic-0.22/trimmomatic-0.22.jar org.usadellab.trimmomatic.TrimmomaticPE -threads 4 -phred33 ./R1_Trimmed.fastq ./R2_Trimmed.fastq R1_paired.fastq R1_unpaired.fastq R2_paired.fastq R2_unpaired.fastq ILLUMINACLIP:/home/clf21/bin/TruSeqAdapters.fa:2:40:15 MINLEN:46

echo Producing QC reports for the raw and processed reads...
fastqc --noextract ./R1.fastq
fastqc --noextract ./R1_paired.fastq

echo Cleaning Up...
gzip *R?_0??.fastq
mv *R?_0??.fastq.gz ./raw/
mv *.csv ./raw/
mv trimlog.txt ./raw/
rm ./R?.fastq
rm ./R1_Trimmed.fastq
rm ./R2_Trimmed.fastq
mkdir tophat || { echo 'ERROR: could not make tophat directory for output, exiting...' ; exit 1; }
mkdir logs || { echo 'ERROR: could not make logs directory, exiting...' ; exit 1; }

echo Beginning Alignment...
tophat -r 150 --mate-std-dev 75 --segment-length 23 --library-type fr-unstranded -p 4 --transcriptome-index=/home/clf21/RNA-seq/GTF_knownGenes/compiled_$Build -x 4 -n 2 -o ./tophat/ /home/clf21/bin/bowtie-0.12.7/indexes/$Build ./R1_paired.fastq ./R2_paired.fastq 2> ./logs/tophat_log.txt
# The -T flag following transcriptome-index has been removed to allow reads to map outside the reference transcriptome
# Might want to change the segment-length parameter based on the extent of read trimming so as not to be greater than double the total read length after trimming

echo Cleaning Up...
rm *paired.fastq
mkdir cufflinks  || { echo 'ERROR: could not make cufflinks directory for output, exiting...' ; exit 1; }

echo Beginning Cufflinks FPKM Estimates...
echo Using Reference Transcriptome found at /home/clf21/RNA-seq/GTF_knownGenes/"$Build"_genes.gtf ...
cufflinks -p 4 -u -o ./cufflinks/ --GTF-guide /home/clf21/RNA-seq/GTF_knownGenes/"$Build"_genes.gtf ./tophat/accepted_hits.bam 2> ./logs/cufflinks_log.txt
# The -G flag for forcing cufflinks estimates on the reference transcriptome was replaced with the GTF-guide flag that will base the estimates off the reference but also assemble new transcripts outside the reference transcriptome where needed.

echo Making bigWig coverage track for browser visualization...
/home/clf21/bin/genomeCoverageBed -split -bg -ibam ./tophat/accepted_hits.bam -g /home/clf21/bin/chrom.sizes."$Build".txt > ./tophat/accepted_hits.bedGraph
bedGraphToBigWig ./tophat/accepted_hits.bedGraph /home/clf21/bin/chrom.sizes."$Build".txt ./tophat/accepted_hits.bigWig
rm ./tophat/accepted_hits.bedGraph
echo Normalizing bigWig signal for comparable browser coverage tracks...
/home/clf21/bin/java_genomics_toolkit/toolRunner.sh wigmath.Scale -i ./tophat/accepted_hits.bigWig -o ./tophat/accepted_hits_norm.wig 2> ./logs/wigScale_log.txt
# By default, this will scale the bigWig coverage to a mean of 1
wigToBigWig ./tophat/accepted_hits_norm.wig /home/clf21/bin/chrom.sizes."$Build".txt ./tophat/accepted_hits_scaled.bigWig
rm ./tophat/accepted_hits_norm.wig
#rm ./tophat/accepted_hits.bigWig


echo 'Done. Run appears to have completed successfully ; Check for any error messages from Tophat/Cufflinks' 

# Create new .txt file (i.e. assemblies.txt) that lists the transcripts.gtf files for all relevant samples
# Run cuffmerge on the reference and assembled transcriptomes like so: cuffmerge -p 4 -g /home/clf21/RNA-seq/GTF_knownGenes/"$Build"_genes.gtf -s /home/clf21/RNA-seq/GTF_knownGenes/Compiled_"$Build".fa  assemblies.txt
# Run cuffdiff on the merged transcriptome and all relevant alignments like so: cuffdiff -o ./cuffdiff/ -b /home/clf21/RNA-seq/GTF_knownGenes/Compiled_"$Build".fa -L Label1,Label2 -u merged_asm/merged.gtf SampleA1/tophat/accepted_hits.bam,SampleA2/tophat/accepted_hits.bam SampleB1/tophat/accepted_hits.bam,SampleB2/tophat/accepted_hits.bam


