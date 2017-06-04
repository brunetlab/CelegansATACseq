#!/bin/bash

if [[ "$#" -lt 2 ]]
then
    echo "USEAGE: Check for enrichment of one set of peaks with respect to a list of other regions in a provided file (1 per line):" 1>&2
    echo "This program accept any bed, or similarly formatted, file and outputs single base Bedtools Fischer enrichments" 1>&2
    echo "$(basename $0) [BED]"  1>&2
    echo "$(basename $1) REGIONS"  1>&2
    echo "$(basename $2) OPTIONAL: OUTPUT_DIR"  1>&2
    echo "   [BED]: File to annotate" 1>&2
    echo "   [REGIONS]: File with one line per region to check for enrichment" 1>&2
    echo "   [OUTPUT_DIR]: Optional" 1>&2
    exit 1
fi

# =============================
# Settings
# =============================

OVERLAP_PORTION="0.5" # that is the portion of a peak that must be intersected by a annotation to be annotated as such

BED=$(echo $1 | sed 's:/$::g')
REGIONS=$(echo $2 | sed 's:/$::g')
if [[ "$#" -gt 2 ]]
then
    cd $3
fi

PEAK_SUFFIX=$(echo "${BED}" | perl -F'\.' -ane 'print $F[-1];')
OFPREFIX=$(basename "${BED}" | sed s/\."${PEAK_SUFFIX}"//g)

# Set up out put
SUMMARY_FILE="${OFPREFIX}_genicFisherEnrichSummary.txt"

echo "# Checking for genic enrichments with ${BED}" > ${SUMMARY_FILE}

while read line
do
	echo ${line} $(sortBed -i ${line} |\
	grep -v 'chrM' - |\
	bedtools fisher -a ${BED} -b - -g /Users/acd13/Desktop/ce10_chromSizes |\
	tail -n1) >> ${SUMMARY_FILE}
done < ${REGIONS}
