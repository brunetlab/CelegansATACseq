#!/usr/bin/python

# Author : Aaron C Daugherty, Stanford University
# 3 Novemeber, 2014

# This script takes 2 bed files of named features and replaces the features in the first file with the features from the second file, as fetermiend by name.
# If the first file contains features not found in th 2nd file, the original feature is reported


from __future__ import print_function
import sys
import os
from optparse import OptionParser
from collections import defaultdict


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

usage = "usage: %prog [options] [input] This will merge named features in a bed file"

opts = OptionParser(usage=usage)

opts.add_option("-i", help="A 6 column bed file to that you want to replace certain features from")

opts.add_option("-r", help="A 6 column bed file that contains the new features to replace with")

options, arguments = opts.parse_args()


# return usage information if no argvs given

if len(sys.argv)<2:

    os.system(sys.argv[0]+" --help")

    sys.exit()



##### INPUTS AND OUTPUTS #####

# name input and outputs

inputfile = open(options.i, 'r')
replaceFile = open(options.r, 'r')


#=================================================
# Go through the replace file and read it into memory
#=================================================
# This dictionary of lists will hold all the information we need
# The outer key is the feature name
# The value is the original line, this allow for multiple entries for a single name

features = defaultdict(list)

while(True):
    line=replaceFile.readline().rstrip()
    info=line.split("\t")
    
    if (len(info) < 6):
	#warning("Done reading ", replaceFile)
	break
    
    features[info[3]].append(line)

replaceFile.close()


#=================================================
# Print out everything
#=================================================
while(True):
    line=inputfile.readline().rstrip()
    info=line.split("\t")
    
    if (len(info) < 6):
	#warning("Done reading ", inputfile)
	break
    
    if (info[3] not in features): # If we haven't seen this feature
	printOut("\t".join(info))
    
    else:
	for replaceLine in features[info[3]]:
	    replaceInfo=replaceLine.split("\t")
	    printOut("\t".join(replaceInfo))
	

#=================================================
# Finish up
#=================================================

inputfile.close()
sys.exit()
