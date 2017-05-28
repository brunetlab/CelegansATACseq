#!/bin/bash
if [[ $# -lt 6 ]]
then
    echo "USAGE: $(basename $0) [modeFile] [oDir] [oPrefix] [origLabelMnemomicsFile] [stateReorderFile] [columnOrderFile]" >&2
    echo 'Will reorder states and add labels' >&2
    echo '[modelFile]: Path to chmm model file' >&2
    echo '[oDir]: output Directory (Will create if it does not exist)' >&2
    echo '[oPrefix]: output file prefix' >&2
    echo '[origLabelMnemomicsFile]: file mapping original state order to mnemonics' >&2
    echo '[stateReorderFile]: mapping original state order to new state order' >&2
    echo '[columnOrderFile]: file containing new order for columns' >&2
    exit 1
fi

CHMM="/Users/acd13/Desktop/ChromHMM/ChromHMM.jar"

MODELFILE=$1
if [[ ! -e ${MODELFILE} ]]
then
    echo "Model file ${MODELFILE} does not exist" >&2
    exit 1
fi

ODIR=$2
[[ ! -d ${ODIR} ]] && mkdir -p ${ODIR}

OPREFIX=$3

LABELMAP=$4
if [[ ! -e ${LABELMAP} ]]
then
    echo "Menomics file ${LABELMAP} does not exist" >&2
    exit 1
fi

STATEORDER=$5
if [[ ! -e ${STATEORDER} ]]
then
    echo "State reorder file ${STATEORDER} does not exist" >&2
    exit 1
fi

COLUMNORDER=$6
if [[ ! -e ${COLUMNORDER} ]]
then
    echo "Column order file ${COLUMNORDER} does not exist" >&2
    exit 1
fi

java -mx4096M -jar $CHMM Reorder -f ${COLUMNORDER} -i ${OPREFIX} -m ${LABELMAP} -o ${STATEORDER} ${MODELFILE} ${ODIR}
