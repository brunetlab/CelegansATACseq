#!/bin/bash
if [[ $# -lt 4 ]]
then
    echo "USAGE: $(basename $0) [chmm_segments_file] [annoDir] [annoName]" >&2
    echo 'Computes overlap enrichment for a set of annotations' >&2
    echo '[chmm_segments_file]: path to chromHMM segments BED file' >&2
    echo '[annoDir]: Directory containing annotation ANCHOR *.txt.gz files to compute neighborhood enrichments against' >&2
    echo '[res]: resolution of segmentation in bp' >&2
    echo '[width]: +/- number of bins you want around anchor point' >&2
    echo '[labelMap] (OPTIONAL): mnemonic label map (with state ordering prefix)' >&2
    exit 1
fi

CHMM="/Users/acd13/Desktop/ChromHMM/ChromHMM.jar"

FILE=$1
if [[ ! -e ${FILE} ]]
then
    echo "Segments file ${FILE} does not exist" >&2
    exit 1
fi

ANNODIR=$2
if [[ ! -d ${ANNODIR} ]]
then
    echo "Annotation directory ${ANNODIR} does not exist" >&2
    exit 1
fi

RES=$3

WIDTH=$4

LABELMAP=''
if [[ $# -gt 4 ]]
then
    LABELMAP=$5
    if [[ ! -e ${LABELMAP} ]]
    then
		echo "LabelMap file ${LABELMAP} does not exist" >&2
	exit 1
    fi
fi

PREFIX=$(echo $(basename $FILE) | sed 's/_segments.*$//g')
IDIR=$(dirname $FILE)

# =====================
# NEIGHBORHOOD ENRICHMENTS
# =====================                                                                                                                                                                                                                                                      

for i in $(find ${ANNODIR} -name '*.txt.gz')
do
    TYPE_PREF=$(basename $i | sed 's/\..*$//g'); 
    if [[ -e ${LABELMAP} ]]
    then
		java -mx4000M -jar $CHMM NeighborhoodEnrichment -b ${RES} -l ${WIDTH} -r ${WIDTH} -m ${LABELMAP} -t ${PREFIX}_${TYPE_PREF}_neigh_fc ${FILE} $i ${IDIR}/${PREFIX}_${TYPE_PREF}_neighborhood
    else
		java -mx4000M -jar $CHMM NeighborhoodEnrichment -b ${RES} -l ${WIDTH} -r ${WIDTH} -t ${PREFIX}_${TYPE_PREF}_neigh_fc ${FILE} $i ${IDIR}/${PREFIX}_${TYPE_PREF}_neighborhood
    fi
done
