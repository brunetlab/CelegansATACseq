#!/bin/bash
if [[ $# -lt 7 ]]
then
    echo "USAGE: $(basename $0) [iDir]  [oDir] [outPrefix] [nstates] [memory]" >&2
    echo 'Will run chromHMM with specified settings' >&2
    echo '[iDir]: directory containing binarized files *_binary*' >&2
    echo '[oDir]: output Directory (Will create if it does not exist)' >&2
    echo '[outPrefix]: output file prefix' >&2
    echo '[nstates]: number of states' >&2
    echo '[memory]: max memory to be used (in GB)' >&2
    echo '[genomeSize]: file containing chromosome sizes' >&2
    echo '[organismID]: e.g. hg19, mm9. Used to compute enrichments' >&2
    exit 1
fi

CHMM="/Users/acd13/Desktop/ChromHMM/ChromHMM.jar"

IDIR=$1
if [[ ! -d ${IDIR} ]]
then
    echo "Directory ${IDIR} does not exist" >&2
    exit 1
fi

#FILELIST=$2
#if [[ ! -e ${FILELIST} ]]
#then
  #  echo "File list ${FILELIST} does not exist" >&2
  #  exit 1
#fi

ODIR=$2
[[ ! -d ${ODIR} ]] && mkdir -p ${ODIR}

OPREFIX=$3

NSTATES=$4

MEM_TEMP=$5
MEM="-mx$((MEM_TEMP * 1024))M"

GSIZE=$6

ORGANISM=$7

java ${MEM} -jar ${CHMM} LearnModel -b 100 -i ${OPREFIX} -l ${GSIZE} ${IDIR} ${ODIR} ${NSTATES} ${ORGANISM}
