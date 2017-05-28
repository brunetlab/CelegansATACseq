#!/bin/bash
if [[ $# -lt 4 ]]
then
    echo "USAGE: $(basename $0) [chmm_segments_file] [annoDir] [annoName]" >&2
    echo 'Computes overlap enrichment for a set of annotations' >&2
    echo '[chmm_segments_file]: path to chromHMM segments BED file' >&2
    echo '[annoDir]: Directory containing annotation BED files to compute overlaps against' >&2
    echo '[annoName]: Name for annotation to use in output file names' >&2
    echo '[res]: Resolution of segmentation' >&2    
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

ANNONAME=$3

RES=$4

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
# OVERLAP ENRICHMENTS
# =====================                                                                                                                                                                                                                                                      
if [[ -e ${LABELMAP} ]]
then
    java -mx4000M -jar $CHMM OverlapEnrichment -b ${RES} -m ${LABELMAP} -t ${PREFIX}_${ANNONAME}_fc ${FILE} ${ANNODIR} ${IDIR}/${PREFIX}_${ANNONAME}_overlap
else
    java -mx4000M -jar $CHMM OverlapEnrichment -b ${RES} -t ${PREFIX}_${ANNONAME}_fc ${FILE} ${ANNODIR} ${IDIR}/${PREFIX}_${ANNONAME}_overlap
fi
