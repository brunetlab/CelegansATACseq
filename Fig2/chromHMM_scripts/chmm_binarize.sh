#!/bin/bash
if [[ $# -lt 7 ]]
then
    echo "USAGE: $(basename $0) [metaDataFileList] [strandShift] [oDir] [GsizeFile] [memory] [chipDir] [poissonThreshold] [controlDir(OPTIONAL)]" >&2
    echo 'Will use chromHMM binarization routine' >&2
    echo '[metaDataFileList]: 4 Column file [cell] [mark] [chipFileName] [controlFileName (OPTIONAL)]' >&2
    echo '[strandShift]: number of bases by which you want to shift reads (fraglen/2)' >&2
    echo '[oDir]: Output directory (Will create one if it doesnt exist)' >&2
    echo '[GsizeFile]: genome size file with lengths of each chromosome [chr] [len]' >&2
    echo '[memory]: max memory in GB to be used' >&2
    echo '[chipDir]: directory containing ChIP alignment files (BED files can be gzipped)' >&2
    echo '[poissonThreshold]:  specifies the tail probability of the poisson distribution that the binarization threshold should correspond to. The default value of this parameter is 0.0001' >&2
    echo '[controlDir] (OPTIONAL): directory containing InputDNA control alignment files (BED files can be gzipped)' >&2
    exit 1
fi

CHMM="/Users/acd13/Desktop/ChromHMM/ChromHMM.jar"

METADATA=$1
if [[ ! -e ${METADATA} ]]
then
    echo "File list ${METADATA} does not exist" >&2
    exit 1
fi

STRANDSHIFT=$2

ODIR=$3
[[ ! -d ${ODIR} ]] && mkdir -p ${ODIR}

GSIZE=$4
if [[ ! -e ${GSIZE} ]]
then
    echo "Genome size file ${GSIZE} does not exist" >&2
    exit 1
fi

MEM_TEMP=$5
MEM="-mx$((MEM_TEMP * 1024))M"

CHIPDIR=$6
if [[ ! -d ${CHIPDIR} ]]
then
    echo "Directory ${CHIPDIR} does not exist" >&2
    exit 1
fi

POISTHRESH=$7

CONTROLDIR=''
if [[ $# -gt 7 ]]
then
    CONTROLDIR=$8
    if [[ ! -d ${CONTROLDIR} ]]
    then
	echo "Directory ${CONTROLDIR} does not exist" >&2
	exit 1
    fi
fi

if [[ -d ${CONTROLDIR} ]]
then
    java ${MEM} -jar $CHMM BinarizeBed -b 100 -c ${CONTROLDIR} -n ${STRANDSHIFT} -p ${POISTHRESH} ${GSIZE} ${CHIPDIR} ${METADATA} ${ODIR}
else
    java ${MEM} -jar $CHMM BinarizeBed -b 100 -n ${STRANDSHIFT} -p ${POISTHRESH} ${GSIZE} ${CHIPDIR} ${METADATA} ${ODIR}
fi
