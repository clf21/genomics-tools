#!/usr/bin/env python2.7

# Written by Thomas Konneker
# This script takes in a bed file from standard input and removes
# entries that map in excess to a particular base pair, given a 
# threshold. Presumably this removes a lot of PCR artifacts.
  
import sys
import string
from optparse import OptionParser

## Command Lines arguments

parser = OptionParser()
parser.add_option('-d', '--depth', type='int',
                  dest='depth', default=5,
                  help='threshold for trimming reads that have the \
                        same mapping start position')
parser.add_option('-r', '--reduce', action="store_true", dest="reducebp")

(options, args)=  parser.parse_args()

# set the trim depth 
trimdepth = options.depth
bpcut = options.reducebp

# fasta entry class
class bedEntry:

    def __init__(self, chromosome, start, end, name, score,strand, fullBedLine):
        self.chromosome = str(chromosome)
        self.startposition = int(start)
        self.endposition = int(end)
        self.name = name
        self.score = score
        self.strand = strand
        self.fullBedLine = fullBedLine

# Bed file reader object.  Includes method for parsing and yielding
# lines.
class bedReader:

    def __init__(self, inFile):
        self.inFile = inFile


    def spitLines (self):

        bedLine = None
        
        for line in self.inFile:
            fullBedLine = line.strip('\n\r')
            wholeline = fullBedLine.split('\t')
            if bpcut:
                if wholeline[5] == '+' :
                    trimmedposition = int(wholeline[1]) + 1
                    fullBedLine = wholeline[0] + '\t' + wholeline [1] + '\t' \
                                 + str(trimmedposition) + '\t' + wholeline[3] + '\t' \
                                 + wholeline[4] + '\t' + wholeline[5]
                else: 
                    trimmedposition = int(wholeline[2]) - 1
                    fullBedLine = wholeline[0] + '\t' + str(trimmedposition) + '\t' \
                                 + wholeline[2] + '\t' + wholeline[3] + '\t' \
                                 + wholeline[4] + '\t' + wholeline[5]
            bedLine = bedEntry(wholeline[0], wholeline[1],
                               wholeline[2], wholeline[3],
                               wholeline[4], wholeline[5],
                               fullBedLine)
            yield bedLine

# Method that takes a bedReader and depth and returns only lines
# that are below the depth cutoff for a given genomic coordinate
	
def trimulate(bedIn, depth):
    bedIn, depth = bedIn, depth
    cutoff = 0
    position = 0
    for line in bedIn.spitLines():
        if line.startposition != position:
            cutoff = 0
            position = line.startposition
            print line.fullBedLine
        else:
            cutoff += 1
            if cutoff >= depth: continue
            else: print line.fullBedLine

# Main Method.  Operates on stdin.
if __name__ == '__main__':
    entry = bedReader(sys.stdin)
    trimulate(entry, trimdepth)

	    

