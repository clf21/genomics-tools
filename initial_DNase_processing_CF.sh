#!/bin/bash

if [ -z "$1" ]; then 
        echo usage: $0 '<single_directory_of_raw_reads> <genome_build> <totalreadlength>'
    exit
fi

WD=$1
Build=$2
ReadLength=$3
Trim=$(($ReadLength-20))
echo Genome Build Set to: $Build
echo Processing files in : $WD
echo Read length set to: $ReadLength , will trim off $Trim bp from 3-prime end of read...

cd $WD || { echo ERROR: could not find $WD , exiting... ; exit 1; }
mkdir raw || { echo 'ERROR: could not make raw directory, exiting...' ; exit 1; }
echo Concatenating raw reads...
cat  *R1_0??.fastq.gz > ./all.fastq.gz || { echo 'Reads appear to have suffix other than fastq.gz, trying sequence.txt.gz ...' ; cat *sequence.txt.gz > ./all.fastq.gz || echo 'Trying R1.fastq.gz ...' ; cat *R1.fastq.gz > ./all.fastq.gz || { echo 'Cannot recognize read suffix, exiting!! ' ; exit 1; } }
gunzip all.fastq.gz

echo Cleaning Up...
mv *R?_0??.fastq.gz ./raw/
mv *sequence.txt.gz ./raw/
mv *.csv ./raw/

echo 'Producing QC report for the raw reads...'
fastqc --noextract ./all.fastq
mkdir bowtie
mv ./all_fastqc.zip ./bowtie/

echo Beginning Bowtie Alignment to $Build ...
# Bowtie is being run, trimming the 3-prime end from each read down to 20bp and searching for the single best alignment for each read with a reduced seed match length of 20, allowing no  mismatches and and no multi-reads (unique alignments only).  Output is SAM format which is directly converted to BAM format for space.
(bowtie --trim3 $Trim -l 20 -n 0 -m 1 -S -p 6 --un unaligned1.fastq -q /home/clf21/bin/bowtie/indexes/$Build ./all.fastq | samtools view -bS -o bowtie1.bam - ) 2> bowtie_run1.txt
(bowtie --trim3 $Trim -l 20 -n 1 -m 1 -S -p 6 --un unaligned2.fastq -q /home/clf21/bin/bowtie/indexes/$Build ./unaligned1.fastq | samtools view -bS -o bowtie2.bam - ) 2> bowtie_run2.txt
samtools merge bowtie.bam bowtie1.bam bowtie2.bam
rm unaligned?.fastq
rm bowtie?.bam
echo 'Sorting bam file and extracting aligned reads...'
samtools sort bowtie.bam bowtie_sorted || { echo 'ERROR: could not sort .bam output from bowtie, exiting...' ; exit 1; }
samtools view -F 4 -b -h -o accepted_hits.bam bowtie_sorted.bam || { echo 'ERROR: could not remove unmapped reads from .bam alignment, exiting...' ; exit 1; } # Remove unmapped reads from output
#samtools view -f 4 -b -h -o unmapped_reads.bam bowtie_sorted.bam
bamToBed -i accepted_hits.bam > sequence.bed || { echo 'ERROR: could not produce sequence.bed from .bam alignment, exiting...' ; exit 1; }

echo 'Cleaning up...'
rm all.fastq
rm bowtie.bam
rm bowtie_sorted.bam
mv unmapped_reads.bam ./bowtie/
mv accepted_hits.bam ./bowtie/
mv bowtie_run.txt ./bowtie/
mv sequence.bed ./bowtie/
mv bowtie_run?.txt ./bowtie/

#echo 'Filtering out potential PCR artifacts (greater than 100 reads mapping to exact same location) and trimming reads to single 5-prime bp...'
#cat ./bowtie/sequence.bed | readTrimmerPlus.py -d 100 -r > ./bowtie/sequence.filt1.bed

echo 'Filtering out potential PCR artifacts (greater than 70% of reads within 31 bp window map to same bp position with at least 5 reads)...'
cat ./bowtie/sequence.bed | windowTrimmer_minLimit.py > ./bowtie/sequence.filt.bed 

echo 'Filtering out ENCODE blacklisted regions and reads that map to mitochondrial genome...'
intersectBed -v -a ./bowtie/sequence.filt.bed -b /home/clf21/bin/"$Build"_blacklist_regions_2014ENCODE.bed | grep -v 'chrM' | grep -v 'chrY' > ./bowtie/sequence.final.bed
diff ./bowtie/sequence.bed ./bowtie/sequence.final.bed > ./bowtie/filtered_reads.bed
rm ./bowtie/sequence.bed
rm ./bowtie/sequence.filt.bed

echo 'Running MACS2 to call peaks and generate bedGraph...'
#note this uses parameters suggested by MACS2 author specifically for DNase 5' end data. Shifting the reads by negative 100bp and then extending 3' end to make reads appear 200bp in length total. This centers the reads on the 5' cut site (+ strand)
mkdir macs2 || { echo 'ERROR: could not make macs2 output directory, exiting...' ; exit 1; }
macs2 callpeak -t ./bowtie/sequence.final.bed -f BED -g hs -n macs2 --outdir ./macs2/ -q 0.05 --nomodel --shift -100 --extsize 200 --bdg 2> ./macs2/macs2_run.txt

echo 'Counting mapped and filter-passing reads for normalization...'
NR=$(wc -l ./bowtie/sequence.final.bed | awk {'print $1'})	# obtaining total number of post-filtering reads
NRM=$(echo "scale=6; $NR/1000000" | bc)	# dividing by 1 million reads
SF=$(echo "scale=6; 1/$NRM" | bc)	# taking the inverse as a scaling factor (this value gets multiplied by raw read coverage at each position in wig file ; results in reads/million mapped)

echo 'Converting bedGraph to bigWig for browser visualization...'
#sort -k1,1 -k2,2n ./macs2/macs2_treat_pileup.bdg > ./macs2/macs2_sorted.bedGraph
bedGraphToBigWig ./macs2/macs2_treat_pileup.bdg /home/clf21/bin/chrom.sizes."$Build".txt ./macs2/macs2.bigWig
rm ./macs2/macs2_treat_pileup.bdg
rm ./macs2/macs2_control_lambda.bdg

echo 'Normalizing bigWig to reads per million aligned...'
/home/clf21/bin/java_genomics_toolkit/toolRunner.sh wigmath.Scale -m "$SF" -i ./macs2/macs2.bigWig -o ./macs2/macs2_scaled.wig
wigToBigWig ./macs2/macs2_scaled.wig /home/clf21/bin/chrom.sizes."$Build".txt ./macs2/macs2_scaled.bigWig
rm ./macs2/macs2_scaled.wig

echo 'Cleaning up...'
gzip ./bowtie/sequence.final.bed
gzip ./bowtie/filtered_reads.bed

echo 'Done.  Please check for any error messages.'


####Deprecated code below:###

### This is an alternative method for removing adapter sequences - it works well but the fasta file has to be added to for barcoded DNase adapters (varying lengths) and it takes a while to run.  
#echo Trimming out DNase Adapter Sequences...
#java -classpath /home/clf21/bin/Trimmomatic-0.22/trimmomatic-0.22.jar org.usadellab.trimmomatic.TrimmomaticSE -threads 6 -phred33 ./R1.fastq R1_trimmed.fastq ILLUMINACLIP:/home/clf21/bin/DNaseAdapters.fa:2:40:15 CROP:23 MINLEN:18 
# This requires reads be between 18 and 23 base pairs in length and removes any near-perfect matches to the adapter sequences.

#echo Making exactCut .wig file for counts data...
#/home/clf21/bin/exactcuts_helper_CF_2.sh $WD $Build || { echo 'ERROR: could not produce exactCuts wiggle files, exiting...' ; exit 1; }
# Note: this exactCuts script currently only has hg19 and mm9 available to it.

#echo Running F-seq to call peaks...
#mkdir peaks_initial_"$Build"
#/nfs/furey_sata2/bin/fseq -v -f 0 -o ./peaks_initial_"$Build"/ -of npf ./bowtie/sequence.final.bed

#echo Running F-seq to produce smoothened coverage bigWig...
#mkdir parzen_"$Build"
#/nfs/furey_sata2/bin/fseq -v -f 0 -o ./parzen_"$Build"/ -of wig ./bowtie/sequence.final.bed
# Note: if we want to, we can take into account copy number differences and mappability differences for a particular species and cell type by specifying -b and -p with files located at /nfs/furey_sata2/uniqueness/...
# The -f 0 parameter is set to work with reads that have been trimmed down to single 5-prime bp

#echo Making exactCuts type of bigWig coverage track for browser visualization...
#mkdir exact_cuts_"$Build"
#/home/clf21/bin/genomeCoverageBed -bg -i ./bowtie/sequence.final.bed -g /home/clf21/bin/chrom.sizes."$Build".txt > ./bowtie/sequence.final.bedGraph
#bedGraphToBigWig ./bowtie/sequence.final.bedGraph /home/clf21/bin/chrom.sizes."$Build".txt ./exact_cuts_"$Build"/sequence.final.bigWig
#rm ./bowtie/accepted_hits.bedGraph

#echo Converting Fseq output wigs to bigWig...
#cat ./parzen_"$Build"/*.wig | wigToBigWig stdin -clip /home/clf21/bin/chrom.sizes."$Build".txt ./parzen_"$Build"/fseq_parzen.bigWig

#echo Creating top 100,000 and 50,000 peak files by Fseq score...
#cat ./peaks_initial_"$Build"/chr*.npf | sort -k7 -n -r | head -n 100000 > ./peaks_initial_"$Build"/top100k.npf
#cat ./peaks_initial_"$Build"/chr*.npf | sort -k7 -n -r | head -n 50000 > ./peaks_initial_"$Build"/top50k.npf

#echo Normalizing bigWig signal for comparable browser coverage tracks...
#/home/clf21/bin/java_genomics_toolkit/toolRunner.sh wigmath.Scale -i ./bowtie/accepted_hits.bigWig -o ./bowtie/accepted_hits_norm.wig
#wigToBigWig ./bowtie/accepted_hits_norm.wig /home/clf21/bin/chrom.sizes."$Build".txt ./bowtie/accepted_hits_scaled.bigWig
#rm ./bowtie/accepted_hits_norm.wig
#echo Producing moving average-smoothed coverage for visualization...
#/home/clf21/bin/java_genomics_toolkit/toolRunner.sh wigmath.MovingAverageSmooth -i ./bowtie/accepted_hits_scaled.bigWig -o ./bowtie/accepted_hits_MA.wig
#wigToBigWig ./bowtie/accepted_hits_MA.wig /home/clf21/bin/chrom.sizes."$Build".txt ./bowtie/accepted_hits_MA.bigWig
#rm ./bowtie/accepted_hits_MA.wig



