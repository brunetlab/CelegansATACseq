#!/usr/bin/python



# Author: Aaron Daugherty, Stanford University
# this is a slightly different version from the original python script in that this doesn't take a bed file, but instead a tsv with the peak names in the first column
# this is done because now I'm splitting the peaks into sections, but I want everything to remain with the peak name
# and example input line:
#chrI_3887_4099	gei-11_M2337_1.02	1
# this is to work in collaboration with /Users/acd13/Desktop/ATAC/Analysis/predictingAccessChangeWithMotifs/getMotifCounts_inSplitPeaks.sh (before)
# and /Users/acd13/Desktop/ATAC/Analysis/predictingAccessChangeWithMotifs/predictingAccessibilityChangesWithMotifs_RF_withSplitPeaks.R (after)
# this isn't the most memory efficient, but it's not big enough to matter


##### IMPORT MODULES #####

# import necessary for python

from __future__ import print_function

import os

import re

import sys

import string

from collections import defaultdict

from optparse import OptionParser


#=================================================
# Define functions
#=================================================

def warning(*objs):
    print("WARNING: ", *objs, file=sys.stderr)

def useage(*objs):
    print("USEAGE: ", *objs, file=sys.stderr)
    
def printOut(*objs):
    print(*objs, file=sys.stdout)

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

#### OPTIONS ####

# define options

opts = OptionParser()

usage = "usage: %prog [options] [inputs] T"

opts = OptionParser(usage=usage)

opts.add_option("-f", help="<motif counts File> input file")

options, arguments = opts.parse_args()


# return usage information if no argvs given

if len(sys.argv) < 1:

    os.system(sys.argv[0]+" --help")

    sys.exit()



##### INPUTS AND OUTPUTS #####

# name input and outputs

inputfileName = options.f

inputfile_rds = open(inputfileName,'r')


##### SCRIPT #####
matrix = defaultdict(list)
# the key is the motif, and the value is a list of peaks

while 1:
    
    # read lines
    
    in_line = inputfile_rds.readline()
    
    # break if at end of file
    
    if not in_line:
        break
    
    fields = in_line.rstrip().split("\t")
    
    if fields[1] in matrix.keys(): # the key is the motif name
        matrix[fields[1]].append("_".join((fields[0],fields[2])))
    else:
        matrix[fields[1]]=["_".join((fields[0],fields[2]))]
inputfile_rds.close()

#=================================================
# Finish up
#========================================

# print out a csv of peaks for each motif
for motif in matrix.keys():
    printOut("\t".join((motif, ",".join(matrix[motif]))))


sys.exit()
