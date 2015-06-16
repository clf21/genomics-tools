#!/usr/bin/env python2.7

Usage = """ Usage: extractFromBamorBed_CF.py <bed_file_of_targets> <bam_or_bed_file_with_mapped_reads> <shift_size> <optional=findMeans>
 The bed file of target regions should have chrom, start, stop, and strand columns. Must have bedTools installed."""

import sys
import subprocess
from numpy import *

if len(sys.argv) < 4:
	print Usage
	sys.exit()

targets = str(sys.argv[1])
counts = str(sys.argv[2])
shiftsize = int(sys.argv[3])
findmeans = str(sys.argv[4])

if shiftsize == 0:
	print "shiftsize set to zero, proceeding to extract scores..."
	print "extracting scores from ", counts, " for each interval found in ", targets, "..." 
	endings = [".bed", ".bed.gz"]
	if counts.lower().endswith('.bam'):
		covBedCommand = str('coverageBed -abam ' + counts + ' -b ' + targets + ' -d > out1.bed')
		print covBedCommand
		subprocess.call(covBedCommand, shell=True)

	elif counts.lower().endswith(tuple(endings)):
		covBedCommand = str('coverageBed -a ' + counts + ' -b ' + targets + ' -d > out1.bed')
		print covBedCommand
		subprocess.call(covBedCommand, shell=True)

	else:
		print "bam or bed file with mapped reads was not recognized, exiting!"
		sys.exit()

else:
	print "shiftsize set to ", str(shiftsize), ", now shifting reads to account for insert size..."
	
	if counts.lower().endswith('.bam'):
		makeBedCommand = str('bamToBed -i ' + counts + ' > tmp_counts.bed')
		print makeBedCommand
		subprocess.call(makeBedCommand, shell=True)
                print "shifting reads ..."
                BedFile = open('tmp_counts.bed', 'rU')
		ShiftedBedFile = open('tmp_counts2.bed', 'w')
                for Line in BedFile:
                        Line = Line.strip('\n')
                        Chr, Start, Stop, Name, Score, Strand = Line.split('\t')
                        if Strand == '+':
                                NewStart = int(Start) + int(shiftsize)
                                NewStop = int(Stop) + int(shiftsize)
                        else:
                                NewStart = int(Start) - int(shiftsize)
                                NewStop = int(Stop) - int(shiftsize)
			if int(NewStart) > 0:
                        	ShiftedBedFile.write(Chr + '\t' + str(NewStart) + '\t' + str(NewStop) + '\t' + str(Name) + '\t' + str(Score) + '\t' + str(Strand) + '\n')

                BedFile.close()
		ShiftedBedFile.close()	

	elif counts.lower().endswith('.bed'):
		print ".bed file supplied, skipping conversion to .bed ..."
		print "shifting reads ..."
		BedFile = open(counts, 'rU')
		ShiftedBedFile = open('tmp_counts2.bed', 'w')
		for Line in BedFile:
			Line = Line.strip('\n')
			Chr, Start, Stop, Name, Score, Strand = Line.split('\t')
			if Strand == '+':
				NewStart = int(Start) + int(shiftsize)
				NewStop = int(Stop) + int(shiftsize)
			else:
				NewStart = int(Start) - int(shiftsize)
				NewStop = int(Stop) - int(shiftsize)
			if int(NewStart) > 0:
				ShiftedBedFile.write(Chr + '\t' + str(NewStart) + '\t' + str(NewStop) + '\t' + str(Name) + '\t' + str(Score) + '\t' + str(Strand) + '\n')
		
		BedFile.close()
		ShiftedBedFile.close()

	else:
		print "bam or bed file with mapped reads was not recognized, exiting!"
		sys.exit()

	print "extracting (shifted) read coverage for .bed coordinates..."
	covBedCommand = str('coverageBed -a ' + 'tmp_counts2.bed'  + ' -b ' + targets + ' -d > out1.bed')
        print covBedCommand
        subprocess.call(covBedCommand, shell=True)

print "reformatting coverageBed output..."
grpByCommand = str("groupBy -i out1.bed -g 1,2,3,4 -c 6 -o collapse | sed 's/,/\t/g' > out2.bed")
print grpByCommand
subprocess.call(grpByCommand, shell=True)

def ReverseStrand():
	CountFile = open("out2.bed", 'rU')
	RevCountFile = open("out3.bed", 'w')
	for Line in CountFile:
    		Line=Line.strip('\n')
    		ElementList = Line.split('\t')
        	Chrom = ElementList[0]
        	Start = ElementList[1]
        	Stop = ElementList[2]
        	Strand = ElementList[3]
        	Scores = ElementList[4:-1]
        	Scores.append(ElementList[-1])
        	if Strand=="-":
                	Scores.reverse()
        	Entry = Chrom + '\t' + str(Start) + '\t' + str(Stop) + '\t' + str(Strand) + '\t' +  "\t".join(Scores) + '\n'
        	if "n/a" in Scores:
				pass
        	else:
				RevCountFile.write(Entry)
	CountFile.close()
	RevCountFile.close()

print "reversing any negative strand interval counts..."
ReverseStrand()

def FindColMeans():
	MeansFile = open("means.txt", 'a')
	data = genfromtxt('out3.bed')
	data = data[:,4:]
	means = data.mean(axis=0)
	print('\t'.join(map(str,means)) + '\n')
	MeansFile.write(counts + '\t' + '\t'.join(map(str,means)) + '\n')
	MeansFile.close()

if findmeans == 'findMeans':
	print "aggregate coverage (mean coverage / bp) for your regions written to means.txt and is: "
	FindColMeans()
else:
	renameTmp = str('mv out3.bed read_counts_out.bed')
	subprocess.call(renameTmp, shell=True)
	print "individual read counts per bp for your regions written to read_counts_out.bed."
		
print "Removing intermediate files..."
removeTmp = str('rm out?.bed ; rm tmp_counts*.bed')
subprocess.call(removeTmp, shell=True)

