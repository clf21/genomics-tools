#!/bin/bash

if [ -z "$1" ]; then 
        echo usage: $0 '<single_directory_of_raw_reads> <barcode_1> <barcode_2> <barcode_3> <barcode_4>    *note that barcodes must be given in numerical order and are given as 1, 2, 4, or 6'
    exit
fi

WD=$1

echo Splitting files found in : $WD
cd $WD || { echo "ERROR: could not find $WD , exiting..." ; exit 1; }
echo "Renaming raw read files if needed..." 
for i in *fastq.gz; do mv "$i" "${i/.fastq.gz}".fastq.sequence.txt.gz || echo "No file renaming needed" ; done

# For splitting 2 barcodes:
if [ $# -eq 3 ]
	then
		BarcodeA=$2
		BarcodeB=$3
		echo Pulling out sequences with barcode number $BarcodeA ...
		if [ $BarcodeA == '1' ] 
    		then
         		zcat *.txt.gz | egrep -A 2 -B 1 '^(.{15,25})AGTCGTGATGTTCG' | grep -v "^--$" | gzip > barcode1_fastq.gz || { echo "ERROR: could not split barcode $BarcodeA, exiting..." ; exit 1; }
		elif [ $BarcodeA == '2' ]
    		then
         		zcat *.txt.gz | egrep -A 2 -B 1 '^(.{15,25})AGTACATCGGTTCG' | grep -v "^--$" | gzip > barcode2_fastq.gz || { echo "ERROR: could not split barcode $BarcodeA, exiting..." ; exit 1; }
		elif [ $BarcodeA == '4' ]
    		then
         		zcat *.txt.gz | egrep -A 2 -B 1 '^(.{15,25})AGTTGGTCAGTTCG' | grep -v "^--$" | gzip > barcode4_fastq.gz || { echo "ERROR: could not split barcode $BarcodeA, exiting..." ; exit 1; }
		elif [ $BarcodeA == '6' ]
    		then
         		zcat *.txt.gz | egrep -A 2 -B 1 '^(.{15,25})AGTATTGGCGTTCG' | grep -v "^--$" | gzip > barcode6_fastq.gz || { echo "ERROR: could not split barcode $BarcodeA, exiting..." ; exit 1; }
		else
    		echo 'ERROR: I do not recognize that barcode! - it must be 1, 2, 4, or 6' ; 
    		exit 1;
		fi

		echo Pulling out sequences with barcode number $BarcodeB ...
		if [ $BarcodeB == '1' ] 
    		then
         		zcat *.txt.gz | egrep -A 2 -B 1 '^(.{15,25})AGTCGTGATGTTCG' | grep -v "^--$" | gzip > barcode1_fastq.gz || { echo "ERROR: could not split barcode $BarcodeB, exiting..." ; exit 1; }
		elif [ $BarcodeB == '2' ]
    		then
         		zcat *.txt.gz | egrep -A 2 -B 1 '^(.{15,25})AGTACATCGGTTCG' | grep -v "^--$" | gzip > barcode2_fastq.gz || { echo "ERROR: could not split barcode $BarcodeB, exiting..." ; exit 1; }
		elif [ $BarcodeB == '4' ]
    		then
         		zcat *.txt.gz | egrep -A 2 -B 1 '^(.{15,25})AGTTGGTCAGTTCG' | grep -v "^--$" | gzip > barcode4_fastq.gz || { echo "ERROR: could not split barcode $BarcodeB, exiting..." ; exit 1; }
		elif [ $BarcodeB == '6' ]
    		then
         		zcat *.txt.gz | egrep -A 2 -B 1 '^(.{15,25})AGTATTGGCGTTCG' | grep -v "^--$" | gzip > barcode6_fastq.gz || { echo "ERROR: could not split barcode $BarcodeB, exiting..." ; exit 1; }
		else
    		echo 'ERROR: I do not recognize that barcode! - it must be 1, 2, 4, or 6' ; 
    		exit 1;
		fi

		echo Pulling out leftover reads without matching barcodes $BarcodeA and $BarcodeB ...
		if [[ $BarcodeA == '1' && $BarcodeB == '2' ]]
    		then
         		zcat *.txt.gz | egrep -v '^(.{15,25})AGTCGTGATGTTCG' |  egrep -v '^(.{15,25})AGTACATCGGTTCG' | egrep -A 2 -B 1 '^[NACTG]{48,}+' | grep -v "^--$" | gzip > unsorted_reads.fastq.gz || { echo "ERROR: could not pull out non-sorted reads, exiting..." ; exit 1; }
		elif [[ $BarcodeA == '4' && $BarcodeB == '6' ]]
    		then
         		zcat *.txt.gz | egrep -v '^(.{15,25})AGTTGGTCAGTTCG' |  egrep -v '^(.{15,25})AGTATTGGCGTTCG' | egrep -A 2 -B 1 '^[NACTG]{48,}+' | grep -v "^--$" | gzip > unsorted_reads.fastq.gz || { echo "ERROR: could not pull out non-sorted reads, exiting..." ; exit 1; }
		elif [[ $BarcodeA == '1' && $BarcodeB == '4' ]]
    		then
    	 		zcat *.txt.gz | egrep -v '^(.{15,25})AGTCGTGATGTTCG' |  egrep -v '^(.{15,25})AGTTGGTCAGTTCG' | egrep -A 2 -B 1 '^[NACTG]{48,}+' | grep -v "^--$" | gzip > unsorted_reads.fastq.gz || { echo "ERROR: could not pull out non-sorted reads, exiting..." ; exit 1; }
		elif [[ $BarcodeA == '2' && $BarcodeB == '4' ]]
    		then
         		zcat *.txt.gz | egrep -v '^(.{15,25})AGTACATCGGTTCG' |  egrep -v '^(.{15,25})AGTTGGTCAGTTCG' | egrep -A 2 -B 1 '^[NACTG]{48,}+' | grep -v "^--$" | gzip > unsorted_reads.fastq.gz || { echo "ERROR: could not pull out non-sorted reads, exiting..." ; exit 1; }
		elif [[ $BarcodeA == '1' && $BarcodeB == '6' ]]
			then
         		zcat *.txt.gz | egrep -v '^(.{15,25})AGTCGTGATGTTCG' |  egrep -v '^(.{15,25})AGTATTGGCGTTCG' | egrep -A 2 -B 1 '^[NACTG]{48,}+' | grep -v "^--$" | gzip > unsorted_reads.fastq.gz || { echo "ERROR: could not pull out non-sorted reads, exiting..." ; exit 1; } 
		elif [[ $BarcodeA == '2' && $BarcodeB == '6' ]]
			then
         		zcat *.txt.gz | egrep -v '^(.{15,25})AGTACATCGGTTCG' |  egrep -v '^(.{15,25})AGTATTGGCGTTCG' | egrep -A 2 -B 1 '^[NACTG]{48,}+' | grep -v "^--$" | gzip > unsorted_reads.fastq.gz || { echo "ERROR: could not pull out non-sorted reads, exiting..." ; exit 1; } 
		else
    		echo 'ERROR: Unrecognized barcode pair, unable to produce unsorted reads!  Check to make sure supplied barcode numbers are in numerical order...' ;
    		exit 1;
		fi

		echo 'Done. Barcoded samples appear to have split successfully.'
fi


# For splitting 3 barcodes:
if [ $# -eq 4 ]
	then
		BarcodeA=$2
		BarcodeB=$3
		BarcodeC=$4
		
		echo Pulling out barcoded sequences $BarcodeA, $BarcodeB, and $BarcodeC ...
		if [[ $BarcodeA == '1' && $BarcodeB == '2' && $BarcodeC == '4' ]]
			then
				zcat *.txt.gz | egrep -A 2 -B 1 '^(.{15,25})AGTCGTGATGTTCG' | grep -v "^--$" | gzip > barcode1_fastq.gz || { echo "ERROR: could not split barcode $BarcodeA, exiting..." ; exit 1; }
				zcat *.txt.gz | egrep -A 2 -B 1 '^(.{15,25})AGTACATCGGTTCG' | grep -v "^--$" | gzip > barcode2_fastq.gz || { echo "ERROR: could not split barcode $BarcodeB, exiting..." ; exit 1; }
				zcat *.txt.gz | egrep -A 2 -B 1 '^(.{15,25})AGTTGGTCAGTTCG' | grep -v "^--$" | gzip > barcode4_fastq.gz || { echo "ERROR: could not split barcode $BarcodeC, exiting..." ; exit 1; }
		elif [[ $BarcodeA == '2' && $BarcodeB == '4' && $BarcodeC == '6' ]]
			then
				zcat *.txt.gz | egrep -A 2 -B 1 '^(.{15,25})AGTACATCGGTTCG' | grep -v "^--$" | gzip > barcode2_fastq.gz || { echo "ERROR: could not split barcode $BarcodeA, exiting..." ; exit 1; }
				zcat *.txt.gz | egrep -A 2 -B 1 '^(.{15,25})AGTTGGTCAGTTCG' | grep -v "^--$" | gzip > barcode4_fastq.gz || { echo "ERROR: could not split barcode $BarcodeB, exiting..." ; exit 1; }
				zcat *.txt.gz | egrep -A 2 -B 1 '^(.{15,25})AGTATTGGCGTTCG' | grep -v "^--$" | gzip > barcode6_fastq.gz || { echo "ERROR: could not split barcode $BarcodeC, exiting..." ; exit 1; }
		elif [[ $BarcodeA == '1' && $BarcodeB == '4' && $BarcodeC == '6' ]]
			then
				zcat *.txt.gz | egrep -A 2 -B 1 '^(.{15,25})AGTCGTGATGTTCG' | grep -v "^--$" | gzip > barcode1_fastq.gz || { echo "ERROR: could not split barcode $BarcodeA, exiting..." ; exit 1; }
				zcat *.txt.gz | egrep -A 2 -B 1 '^(.{15,25})AGTTGGTCAGTTCG' | grep -v "^--$" | gzip > barcode4_fastq.gz || { echo "ERROR: could not split barcode $BarcodeB, exiting..." ; exit 1; }
				zcat *.txt.gz | egrep -A 2 -B 1 '^(.{15,25})AGTATTGGCGTTCG' | grep -v "^--$" | gzip > barcode6_fastq.gz || { echo "ERROR: could not split barcode $BarcodeC, exiting..." ; exit 1; }
		else
    		echo 'ERROR: I do not recognize those barcodes! - each must be 1, 2, 4, or 6' ; 
    		exit 1;
    	fi
    	
    	echo Pulling out leftover reads without matching barcodes 1, 2, 4, or 6 ...
		zcat *.txt.gz | egrep -v '^(.{15,25})AGTCGTGATGTTCG' |  egrep -v '^(.{15,25})AGTACATCGGTTCG' | egrep -v '^(.{15,25})AGTTGGTCAGTTCG' |  egrep -v '^(.{15,25})AGTATTGGCGTTCG' | egrep -A 2 -B 1 '^[NACTG]{48,}+' | grep -v "^--$" | gzip > unsorted_reads.fastq.gz || { echo "ERROR: could not pull out non-sorted reads, exiting..." ; exit 1; }
		
		echo 'Done. Barcoded samples appear to have split successfully.'
fi

			
# For splitting all 4 barcodes:
if [ $# -eq 5 ]
	then
		BarcodeA=$2
		BarcodeB=$3
		BarcodeC=$4
		BarcodeD=$5
		
		echo Pulling out barcoded sequences 1, 2, 4, and 6 ...
		zcat *.txt.gz | egrep -A 2 -B 1 '^(.{15,25})AGTCGTGATGTTCG' | grep -v "^--$" | gzip > barcode1_fastq.gz || { echo "ERROR: could not split barcode $BarcodeA, exiting..." ; exit 1; }
		zcat *.txt.gz | egrep -A 2 -B 1 '^(.{15,25})AGTACATCGGTTCG' | grep -v "^--$" | gzip > barcode2_fastq.gz || { echo "ERROR: could not split barcode $BarcodeB, exiting..." ; exit 1; }
		zcat *.txt.gz | egrep -A 2 -B 1 '^(.{15,25})AGTTGGTCAGTTCG' | grep -v "^--$" | gzip > barcode4_fastq.gz || { echo "ERROR: could not split barcode $BarcodeC, exiting..." ; exit 1; }
		zcat *.txt.gz | egrep -A 2 -B 1 '^(.{15,25})AGTATTGGCGTTCG' | grep -v "^--$" | gzip > barcode6_fastq.gz || { echo "ERROR: could not split barcode $BarcodeD, exiting..." ; exit 1; }
		
		echo Pulling out leftover reads without matching barcodes 1, 2, 4, or 6 ...
		zcat *.txt.gz | egrep -v '^(.{15,25})AGTCGTGATGTTCG' |  egrep -v '^(.{15,25})AGTACATCGGTTCG' | egrep -v '^(.{15,25})AGTTGGTCAGTTCG' |  egrep -v '^(.{15,25})AGTATTGGCGTTCG' | egrep -A 2 -B 1 '^[NACTG]{48,}+' | grep -v "^--$" | gzip > unsorted_reads.fastq.gz || { echo "ERROR: could not pull out non-sorted reads, exiting..." ; exit 1; }
		
		echo 'Done. Barcoded samples appear to have split successfully.'
fi

