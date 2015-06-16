#!/usr/bin/env python2.7

# Written by Thomas Konneker
# This script is intended to crawl through bed files to filter
# out reads that are too concentrated within a specified base-pair 
# window.  This is presumably because they are PCR artifacts. 
# See Boyle et. al 2011, Genome Research
# Currently this program does not 'window' in the traditional sense
# of tiling out the genome and then evaluating those
# windows. Instead starts a window at the next place there is an entry
# that follows the previous window.

# Modified 7-1-2014 to include requirement for 5 reads minimum for removal of putative artifact.

import sys
import string
from optparse import OptionParser
from collections import Counter

##  Command Line

parser = OptionParser()
parser.add_option('-w', '--windowsize', type='int',
          dest='windowsize', default=31,
          help='threshold for trimming reads that have \
          the same mapping start position')
parser.add_option('-c', '--cutoff', type="float",
          dest='cutoff', default=0.70,
          help="threshold for concentration at a single base \
          within a window to cutoff")

(options, args)=  parser.parse_args()


# set the window size, cutoff from the provided arguments 

windowsize = options.windowsize
cutoff = options.cutoff


# fasta entry class

class bedEntry:

    def __init__(self, chromosome, start, end, name,
                 score, strand, fullBedLine):
        self.chromosome = str(chromosome)
        self.startposition = int(start)
        self.endposition = int(end)
        self.name = name
        self.score = score
        self.strand = strand
        self.fullBedLine = fullBedLine

# bed file reader object.  
# For parsing and yielding bed entries from a file.

class bedReader:

    def __init__(self, inFile):
        self.inFile = inFile

# Method for returning a parsed line of .bed
    
    def spitLines (self):

        bedLine = None
        
        for line in self.inFile:
            fullBedLine = line.strip('\n\r')
            wholeline = fullBedLine.split('\t')
        
            bedLine = bedEntry(wholeline[0], wholeline[1], wholeline [2],
             wholeline[3], wholeline[4], wholeline[5], fullBedLine)
            yield bedLine

# Method to evaluate whether a list of bed is too concentrated at a
# single position, given a concentration threshold. Returns Boolean.

def windowChooser(position_list, cutoff):
    position_list, cutoff = position_list, cutoff
    if position_list:
        count = Counter(position_list)
        peak = float(count.most_common()[0][1])
        if (peak/len(position_list) >= cutoff and peak >= 5): return 0
        else: return 1
    else: return 1


# This method takes a bed file, window size, and concentration variables
# and returns the bedlines that fall under the cutoff criteria for 
# concentration of reads to a single base pair within a window.
    
def windowizer(bedIn, windowsize, cutoff):
    bedIn, windowsize, cutoff = bedIn, windowsize, cutoff
    
    windowstart = 0
    windowend = 0
    current_chr = "chr1"
    full_bed_list = []
    position_list = []
    for line in bedIn.spitLines():
        if line.startposition > windowend:
            windowstart = line.startposition
            windowend = windowstart + windowsize
            if windowChooser(position_list, cutoff):
                for item in full_bed_list:print item
            position_list = [line.startposition]
            full_bed_list = [line.fullBedLine]    
        elif line.startposition <= windowend:
            if line.chromosome == current_chr:
                position_list.append(line.startposition)
                full_bed_list.append(line.fullBedLine)
            else:
                current_chr = line.chromosome
                windowstart = line.startposition
                windowend = windowstart + windowsize
                if windowChooser(position_list, cutoff):
                    for item in full_bed_list:print item
                position_list = [line.startposition]
                full_bed_list = [line.fullBedLine]

        else:
            if windowChooser(position_list, cutoff):
                for item in full_bed_list: print item

# Main Method. Takes from standard in.

if __name__ == '__main__':
    entry = bedReader(sys.stdin)
    windowizer(entry, windowsize, cutoff)
        
