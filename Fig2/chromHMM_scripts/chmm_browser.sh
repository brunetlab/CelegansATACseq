#!/bin/bash
if [[ $# -lt 5 ]]
then
    echo "USAGE: $(basename $0) [segmentFile] [oDir] [colorMap] [labelMap] [memory]" >&2
    echo 'Will create browser tracks for [segmentFile] in [oDir]' >&2
    echo '[segmentFile]: Path to chmm segment file' >&2
    echo '[oDir]: output Directory (Will create if it does not exist)' >&2
    echo '[colorMap]: colormap file' >&2
    echo '[labelMap]: mneomic label map file' >&2
    echo '[memory]: max memory to be used' >&2
    exit 1
fi

CHMM="/Users/acd13/Desktop/ChromHMM/ChromHMM.jar"

SEGFILE=$1
if [[ ! -e ${SEGFILE} ]]
then
    echo "Segment file ${SEGFILE} does not exist" >&2
    exit 1
fi

ODIR=$2
[[ ! -d ${ODIR} ]] && mkdir -p ${ODIR}

COLORMAP=$3
if [[ ! -e ${COLORMAP} ]]
then
    echo "Color map file ${COLORMAP} does not exist" >&2
    exit 1
fi

LABELMAP=$4
if [[ ! -e ${LABELMAP} ]]
then
    echo "Color map file ${LABELMAP} does not exist" >&2
    exit 1
fi

MEM_TEMP=$5
MEM="-mx$((MEM_TEMP * 1024))M"

SEGNAME=$( echo $(basename ${SEGFILE}) | sed -r 's/_segments\.bed//g');
echo $SEGNAME

OPREFIX="${ODIR}/${SEGNAME}"

java ${MEM} -jar $CHMM MakeBrowserFiles -c ${COLORMAP} -m ${LABELMAP} ${SEGFILE} ${SEGNAME} ${OPREFIX}
