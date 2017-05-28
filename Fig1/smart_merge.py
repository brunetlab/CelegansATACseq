# =======================================
# Smart merge for ATAC-seq data
# Author: Daniel Kim
# 11/6/2014
# =======================================

'''
To use this script, do the following:
1) Set up an summary file of all your input files that you want merged. Each peak file should be on a 
different line, and the summit file should be listed next to it (tab separated). Ie, the file should 
look like:

timepoint1.bed timepoint1.summits
timepoint2.bed timepoint2.summits
etc

Remember to include the complete path to the file for each file in the summary file.

2) Run the script as follows:

python smart_merge.py <summary file> <summit dist> <output file>

where:

<summary file> is the file created in step 1
<summit dist> is the distance between peaks desired (suggest: 300)
<output file> name of your output bed file
'''

# Import statements
import sys
import os

# Read in arguments from command line
inputs = sys.argv[1]
summit_dist = int(sys.argv[2])
output_file = sys.argv[3]

# Pull out the prefix of the file from the inputs
prefix = os.path.splitext(inputs)[0]

# Set up lists:
summit_files = []
peak_files = []

# ===================================================
# Read all files to be merged into the same bed file
# And also add in summit info with peak region
# ===================================================

# Open summary file and an output file
IN = open(inputs, 'r')
OUT = open(prefix+"_temp.bed", 'w')

# Set up the lists of peak files and summit files from the input file
for line in IN:
    mark_files = line.rstrip().split()

    # Open the peak file and summit file simultaneously
    PEAK = open(mark_files[0], 'r')
    peak_files.append(PEAK)

    SUMMIT = open(mark_files[1], 'r')
    summit_files.append(SUMMIT)

    # While there are lines to be read in the files,
    # grab the peak region and summit region and put
    # into same file
    while (True):
        peak_line = PEAK.readline()
        summit_line = SUMMIT.readline()
        
        if peak_line == "":
            break

        # Split line
        peak_info = peak_line.rstrip().split()
        summit_info = summit_line.rstrip().split()

        # Get relevant parts
        chrom, peak_start, peak_stop = peak_info[0], peak_info[1], peak_info[2]
        summit_start, summit_stop = summit_info[1], summit_info[2]

        # Write into output file
        OUT.write(chrom+'\t'+summit_start+'\t'+summit_stop+'\t'+peak_start+'\t'+peak_stop+'\n')


IN.close()
OUT.close()

# ===============================================
# Now sort the giant file using UNIX sort
# ===============================================

command = 'sort -k 1,1 -k2,2n {0} > {1}'.format(prefix+"_temp.bed", prefix+"_temp_sorted.bed")
os.system(command)

# ================================================================
# Now merge by lines: if the peaks are within summit_dist of
# each other, merge the peak regions of those summits. 
# ================================================================

IN = open(prefix+"_temp_sorted.bed", 'r')
OUT = open(output_file, 'w')

start = True
peak_num = 0

# For each line in the file
for line in IN:

    # If you're at the start of the file, read in the first peak
    if start:
        summit_info = line.rstrip().split()
        current_chrom, current_summit, current_start, current_stop = summit_info[0], int(summit_info[1]), int(summit_info[3]), int(summit_info[4])
        start = False

    # Otherwise, compare next peak to current peak
    else:
        summit_info = line.rstrip().split()
        chrom, summit_start, summit_stop, peak_start, peak_stop = summit_info[0], int(summit_info[1]), int(summit_info[2]), int(summit_info[3]), int(summit_info[4])

        # If the peak is within summit_dist of the current peak, then join regions
        if abs(summit_start - current_summit) < summit_dist and chrom == current_chrom:
            current_summit = summit_start
            if peak_start < current_start:
                current_start = peak_start
            if peak_stop > current_stop:
                current_stop = peak_stop
                
        # Otherwise, "close" the region and write out to file
        else:
            OUT.write(current_chrom+'\t'+str(current_start)+'\t'+str(current_stop)+'\tpeak'+str(peak_num)+'\n')
            peak_num = peak_num + 1
            current_chrom, current_summit, current_start, current_stop = summit_info[0], int(summit_info[1]), int(summit_info[3]), int(summit_info[4])
            
IN.close()
OUT.close()

# Remove temporary files
command='rm {0}_temp*'.format(prefix)
os.system(command)
