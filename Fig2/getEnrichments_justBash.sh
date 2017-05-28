#!/bin/bash

if [[ "$#" -lt 4 ]]
then
    echo "USAGE: $(basename $0) [PEAKS] [SHUFFLED_PEAKS] [LIST_OF_REGIONS] [OVERLAP]" 1>&2
    echo 'This will use a custom perl script to count the number of intersections between a bed file of peaks, and any number of regions of interest, as well as the null background for these intersections' 1>&2
    echo "   [PEAKS]: Bed file of peaks" 1>&2
    echo "   [SHUFFLED_PEAKS_DIR]: The full path of the directory with  shuffled beds, should be gzipped and can be created with createShuffledBedDirectory.sh" 1>&2
    echo "   [LIST_OF_REGIONS]: A text file listing the full path to the regions of interest (also bed files); 1 line per region" 1>&2
    echo "   [OVERLAP]: The portion of the peaks that must overlap to be considered a hit (>0 OVERLAP <=1)" 1>&2
    exit 1
fi

PEAKS=$(echo $1 | sed 's:/$::g')

PEAK_PREFIX=$(basename "${PEAKS}")

SHUFFLED_PEAKS_DIR=$(echo $2 | sed 's:/$::g')

LIST_OF_REGIONS=$(echo $3 | sed 's:/$::g')

OVERLAP=$(echo $4 | sed 's:/$::g')

# Now go through the regions and actually get the enrichments

function anaylzeShuffled {
	REGION_FILE=$1
	shuffledPeaks=$2
	overlap=$3
	outFile=$4
	recipOutFile=$5
	intersectBed -a $shuffledPeaks -b $REGION_FILE -u -f $overlap | wc -l | tr -d ' ' >> $outFile
	intersectBed -a $REGION_FILE -b $shuffledPeaks -u -f $overlap | wc -l | tr -d ' ' >> $recipOutFile
}
export -f anaylzeShuffled

while read line ;do
	# Get the next region
	REGION=$(echo $line)
	echo $REGION
	# Make sure it's actually a file
	if [[ ! -f ${REGION} ]]
	then
            echo "Region file ${REGION} does not exist" >&2
    	    exit 1
	fi

	REGION_PREFIX=$(basename "${REGION}")
	outFile="${REGION_PREFIX}_peaks_${OVERLAP}intersecting_${PEAK_PREFIX}.txt"
	recipOutFile="${REGION_PREFIX}_peaks_${OVERLAP}intersecting_${PEAK_PREFIX}_recip.txt"

	# First count the total number of features we're working with here
	peakNum=$(zcat -f ${PEAKS} | wc -l | tr -d ' ')
	regionNum=$(zcat -f ${REGION} | wc -l | tr -d ' ')
	# Find the number of peaks that intersect the region file
	peakInRegion=$(intersectBed -a $PEAKS -b $REGION -u -f $OVERLAP | wc -l | tr -d ' ')
	regionInPeak=$(intersectBed -a $REGION -b $PEAKS -u -f $OVERLAP | wc -l | tr -d ' ' )
	echo "#${peakNum}" > ${outFile}
	echo "#${peakInRegion}" >> ${outFile}
	echo "#${regionNum}" > ${recipOutFile}
	echo "#${regionInPeak}" >> ${recipOutFile}

	# this is dumb, but parallel won't take 10,000 things, so we split it up
	for i in {1..9}; do
		parallel 'anaylzeShuffled {}' ::: ${REGION} ::: ${SHUFFLED_PEAKS_DIR}/*shuffled.${i}*gz ::: ${OVERLAP} ::: ${outFile} ::: ${recipOutFile}
	done

done < ${LIST_OF_REGIONS}
