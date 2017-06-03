#!/usr/bin/python

# To do:
# Currently the score of the first feature is reported, and any features merged into that feature do not have their score reported

# Author : Aaron C Daugherty, Stanford University
# 13 November, 2015

# This script take a fixed step wig file and converts it to a bedgraph, keeping the step size as the interval

from __future__ import print_function
import sys
import os
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


#=================================================
# Read in the arguments/define options
#=================================================

opts = OptionParser()

usage = "usage: %prog [options] [input] This will convert wigs to bedgraphs, but note that the last step at the 3' end of each chrosome will spill over"

opts = OptionParser(usage=usage)

opts.add_option("-w", help="A fixed step wig file")

options, arguments = opts.parse_args()


# return usage information if no argvs given

if len(sys.argv)<1:

    os.system(sys.argv[0]+" --help")

    sys.exit()



##### INPUTS AND OUTPUTS #####

# name input and outputs

inputfile = open(options.w, 'r')

#=================================================
# Go through the file and print it out directly
#=================================================

chrom = start = step = span = None

while(True):
    info=inputfile.readline().rstrip().split()
    
    if not info:
        break
    
    # we've found a new chromosome, save the necessary data
    if info[0] == "fixedStep":
        chrom = info[1].split("=")[1]
        start = int(info[2].split("=")[1])
        step = int(info[3].split("=")[1])
        span = int(info[4].split("=")[1])
    else:
        if chrom is None:
            warning("First line does not start with fixedStep, which is the format required")
            sys.exit()
        printOut("\t".join([chrom, str(start - 1), str(start - 1 + span), info[0]]))
        start = start + step



#=================================================
# Finish up
#=================================================

inputfile.close()
sys.exit()
