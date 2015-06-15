#!/bin/bash

if [ -z "$1" ]; then 
        echo usage: $0 '<single_directory_of_raw_reads> <genome_build> '
    exit
fi

WD=$1
Build=$2
echo Genome Build Set to: $Build
echo Processing files in : $WD

cd $WD || { echo ERROR: could not find $WD , exiting... ; exit 1; }

if [ -d "./bowtie" ]; then
	echo 'ERROR: Looks like this sample has already been processed, exiting...'; exit 1
fi

echo 'Moving old pipeline-processed data to old pipeline directory...'
mkdir old_pipeline || { echo ERROR: could not write new directory old pipeline, exiting ... ; exit 1; }
mv bowtie ./old_pipeline/
mv macs ./old_pipeline/
mv QC ./old_pipeline/

mv raw/*.fastq.gz ./
cat *.fastq.gz > all.fastq.gz || { echo 'ERROR: did not find fastq files in supplied directory, exiting...' ; exit 1; }
#mkdir raw || { echo 'ERROR: could not make raw directory, exiting...' ; exit 1; }
gunzip all.fastq.gz || { echo 'ERROR: could not unzip fastq file for alignment, exiting...' ; exit 1; }
mv *.fastq.gz ./raw/
mv *.csv ./raw/

echo 'Producing QC report for the raw reads...'
/home/clf21/bin/fastqc --noextract ./all.fastq
mkdir bowtie
mv ./all_fastqc.zip ./bowtie/

echo Beginning Bowtie Alignment to $Build ...
# The first alignment outputs reads that align uniquely to one genomic region with zero mismatches in the seed - usually this is around 60-70% of total reads input
(/home/clf21/bin/bowtie/bowtie -n 0 -m 1 -S -p 6 --un unaligned1.fastq -q /home/clf21/bin/bowtie/indexes/$Build ./all.fastq | samtools view -bS -o bowtie1.bam - ) 2> bowtie_run1.txt
# The second alignment outputs reads that align uniquely to one genomic region with 1 mismatch (such as a SNP) in the seed - this is usually an additional 1-5% of total reads
(/home/clf21/bin/bowtie/bowtie -n 1 -m 1 -S -p 6 --un unaligned2.fastq -q /home/clf21/bin/bowtie/indexes/$Build ./unaligned1.fastq | samtools view -bS -o bowtie2.bam - ) 2> bowtie_run2.txt
samtools merge bowtie.bam bowtie1.bam bowtie2.bam
rm unaligned?.fastq # comment this out if you want to look at unaligned reads
rm bowtie?.bam
echo 'Sorting bam file and extracting aligned reads...'
samtools sort bowtie.bam bowtie_sorted || { echo 'ERROR: could not sort .bam output from bowtie, exiting...' ; exit 1; }
samtools view -F 4 -b -h -o accepted_hits.bam bowtie_sorted.bam || { echo 'ERROR: could not remove unmapped reads from .bam alignment, exiting...' ; exit 1; }
bamToBed -i accepted_hits.bam > sequence.bed || { echo 'ERROR: could not produce sequence.bed from .bam alignment, exiting...' ; exit 1; }

echo 'Determining fragment length shift and QC statistics with phantompeakqualtools...'
mkdir QC
Rscript /home/clf21/bin/phantompeakqualtools/run_spp.R -c=accepted_hits.bam -savp -out=SPP_QC_stats.txt
mv ./SPP_QC_stats.txt ./QC/
mv accepted_hits.pdf ./QC/

echo 'Cleaning up...'
rm all.fastq
rm bowtie.bam
rm bowtie_sorted.bam
mv accepted_hits.bam ./bowtie/
mv bowtie_run?.txt ./bowtie/
mv sequence.bed ./bowtie/

# This is an alternative filter to remove potential PCR artifacts but this may remove real and useful data and imposes a hard limit on the number of reads mapping to each position - as sequencing depth changes, this hard limit becomes troublesome.  Avoiding this for now, and MACS2 will automatically deal with PCR artifact filtering for peak calling downstream.
#echo 'Filtering out potential PCR artifacts (greater than 20 reads mapping to exact same location)...'
#cat ./bowtie/sequence.bed | /home/clf21/bin/readTrimmerPlus.py -d 20 > ./bowtie/sequence.filt1.bed

# This simple filter is very conservative but removes anything that looks like a random pileup of reads at the same position.
echo 'Filtering out potential PCR artifacts (greater than 70% of reads within 31 bp window map to same bp position and it has more than 5 reads)...'
cat ./bowtie/sequence.bed | /home/clf21/bin/windowTrimmer_minLimit.py > ./bowtie/sequence.filt.bed 

# ENCODE has compiled a list of regions to avoid mapping to for both human and mouse genomes that are known repetitive elements, etc. Also removing any reads that mapped to mitochondrial DNA or chromosome Y here.
echo 'Filtering out ENCODE blacklisted regions and reads that map to mitochondrial genome...'
intersectBed -v -a ./bowtie/sequence.filt.bed -b /home/clf21/bin/"$Build"_blacklist_regions_2014ENCODE.bed | grep -v 'chrM' | grep -v 'chrY' > ./bowtie/sequence.final.bed
diff ./bowtie/sequence.bed ./bowtie/sequence.final.bed > ./bowtie/filtered_reads.bed
rm ./bowtie/sequence.bed
rm ./bowtie/sequence.filt.bed

echo 'Running MACS2 to call initial peaks and generate bedGraph...'
mkdir macs2 || { echo 'ERROR: could not make macs2 output directory, exiting...' ; exit 1; }
if [ "$Build" = "hg19" ]
then 
  /home/ls80/local/bin/macs2 callpeak -t ./bowtie/sequence.final.bed -f BED -g hs -n macs2 --outdir ./macs2/ -q 0.05 --bdg --verbose 3 2> ./macs2/macs2_run.txt
else
  /home/ls80/local/bin/macs2 callpeak -t ./bowtie/sequence.final.bed -f BED -g mm -n macs2 --outdir ./macs2/ -q 0.05 --bdg --verbose 3 2> ./macs2/macs2_run.txt 
fi

echo 'Counting mapped and filter-passing reads for normalization...'
NR=$(wc -l ./bowtie/sequence.final.bed | awk {'print $1'})	# obtaining total number of post-filtering reads
NRM=$(echo "scale=6; $NR/1000000" | bc)	# dividing by 1 million reads
SF=$(echo "scale=6; 1/$NRM" | bc)	# taking the inverse as a scaling factor (this value gets multiplied by raw read coverage at each position in wig file ; results in reads/million mapped)

echo 'Converting bedGraph to bigWig for browser visualization...'
sort -k1,1 -k2,2n ./macs2/macs2_treat_pileup.bdg > ./macs2/macs2_sorted.bdg
slopBed -i ./macs2/macs2_sorted.bdg -g /home/clf21/bin/chrom.sizes."$Build".txt -b 0 | /home/clf21/bin/bedClip stdin /home/clf21/bin/chrom.sizes."$Build".txt ./macs2/macs2_sorted_clipped.bdg
/home/clf21/bin/bedGraphToBigWig ./macs2/macs2_sorted_clipped.bdg /home/clf21/bin/chrom.sizes."$Build".txt ./macs2/macs2.bigWig
rm ./macs2/macs2*.bdg

echo 'Normalizing bigWig to reads per million aligned...'
/home/clf21/bin/java_genomics_toolkit/toolRunner.sh wigmath.Scale -m "$SF" -i ./macs2/macs2.bigWig -o ./macs2/macs2_scaled.wig
/home/clf21/bin/wigToBigWig ./macs2/macs2_scaled.wig /home/clf21/bin/chrom.sizes."$Build".txt ./macs2/macs2_scaled.bigWig
rm ./macs2/macs2_scaled.wig

echo 'Cleaning up...'
gzip ./bowtie/sequence.final.bed
gzip ./bowtie/filtered_reads.bed

echo 'Done. Please check for any error messages.'


###Deprecated code below:###

#echo 'Converting .bam alignment to tagAlign format'
#samtools view -F 0x0204 -o - <bamFile> | awk 'BEGIN{OFS="\t"}{if (and($2,16) > 0) {print $3,($4-1),($4-1+length($10)),"N","1000","-"} else {print $3,($4-1),($4-1+length($10)),"N","1000","+"} }' | gzip -c > <gzip_TagAlignFileName>

#echo 'Determining fragment length shift and QC statistics with phantompeakqualtools...'
#mkdir QC
#Rscript /home/clf21/bin/phantompeakqualtools/run_spp.R -c=./bowtie/accepted_hits.bam -savp -out=SPP_QC_stats.txt
#mv ./SPP_QC_stats.txt ./QC/
#mv ./bowtie/accepted_hits.pdf ./QC/

#echo Making bigWig coverage track for browser visualization...
#/home/clf21/bin/genomeCoverageBed -bg -ibam ./bowtie/accepted_hits.bam -g /home/clf21/bin/chrom.sizes."$Build".txt > ./bowtie/accepted_hits.bedGraph
#bedGraphToBigWig ./bowtie/accepted_hits.bedGraph /home/clf21/bin/chrom.sizes."$Build".txt ./bowtie/accepted_hits.bigWig
#rm ./bowtie/accepted_hits.bedGraph
#echo Normalizing bigWig signal for comparable browser coverage tracks...
#/home/clf21/bin/java_genomics_toolkit/toolRunner.sh wigmath.Scale -m "$SF" -i ./bowtie/accepted_hits.bigWig -o ./bowtie/accepted_hits_scaled.wig
#wigToBigWig ./bowtie/accepted_hits_scaled.wig /home/clf21/bin/chrom.sizes."$Build".txt ./bowtie/accepted_hits_scaled.bigWig
#rm ./bowtie/accepted_hits_scaled.wig




