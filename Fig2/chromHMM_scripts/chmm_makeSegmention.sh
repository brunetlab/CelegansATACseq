#!/bin/bash
if [[ $# -lt 6 ]]
then
    echo "USAGE: $(basename $0) [modeFile] [oDir] [oPrefix] [origLabelMnemomicsFile] [stateReorderFile] [columnOrderFile]" >&2
    echo 'Will make the segmentation bed fies given a model (Does not need to re-learn the model)' >&2
    echo '[modelFile]: Path to chmm model file' >&2
    echo '[iDir]: directory containing binarized files *_binary*' >&2
	echo '[oDir]: output Directory (Will create if it does not exist)' >&2
	echo '[oPrefix]: output file prefix' >&2
	echo '[binSize]: size of bins in bp' >&2
    echo '[genomeSize]: file containing chromosome sizes' >&2
    exit 1
fi

CHMM="/Users/acd13/Desktop/ChromHMM/ChromHMM.jar"

MODELFILE=$1
if [[ ! -e ${MODELFILE} ]]
then
    echo "Model file ${MODELFILE} does not exist" >&2
    exit 1
fi

IDIR=$2
if [[ ! -d ${IDIR} ]]
then
    echo "Directory ${IDIR} does not exist" >&2
    exit 1
fi

ODIR=$3
[[ ! -d ${ODIR} ]] && mkdir -p ${ODIR}

OPREFIX=$4

BIN_SIZE=$5

GSIZE=$6

java -mx4096M -jar $CHMM MakeSegmentation -l ${GSIZE} -b ${BIN_SIZE} -i ${OPREFIX} ${MODELFILE} ${IDIR} ${ODIR}

