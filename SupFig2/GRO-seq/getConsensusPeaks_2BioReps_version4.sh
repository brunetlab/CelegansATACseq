#!/bin/bash

#### This is used to get a consensus peak set from a dataset with 2 biological replicates, that were merged and psuedo-replicates created
# It will take all merged peaks that overlap in BOTH pseudoreps and atleast 1 true rep, but aren't in the blacklisted region of ce10



if [[ "$#" -lt 5 ]]
then
	echo "USEAGE: determine consensus peak set:" 1>&2
	echo "This program accepts peak files (e.g. from MACS2) from a dataset with 3 biological replicates that were merged, and then psuedo-rep'd" 1>&2
	echo "This will identify any peaks in the merged peak file that intersect a peak in both psuedoreps and at least 2 true reps"
    echo "$(basename $0) [MERGED_PEAKS]"  1>&2
    echo "$(basename $1) [BIO_REP_A]"  1>&2
    echo "$(basename $2) [BIO_REP_B]"  1>&2
    echo "$(basename $3) [PR_1]"  1>&2
    echo "$(basename $4) [PR_2]"  1>&2
    echo "   [MERGED_PEAKS]: The peaks called from the merged bioReps" 1>&2
    echo "   [BIO_REP_A]: The peaks called from the first biological replicate" 1>&2
    echo "   [BIO_REP_B]: The peaks called from the second biological replicate" 1>&2
    echo "   [PR_1]: The peaks called from the first psuedoreplicate" 1>&2
    echo "   [PR_2]: The peaks called from the first psuedoreplicate" 1>&2
    exit 1
fi

# =============================
# Settings
# =============================

OVERLAP_PORTION="0.5" # that is the portion of a merged peak that must be intersected by a peak in the comparison files (either the reps or the prs)

MERGED_PEAKS=$(echo $1 | sed 's:/$::g')
BIO_REP_A=$(echo $2 | sed 's:/$::g')
BIO_REP_B=$(echo $3 | sed 's:/$::g')
PR_1=$(echo $4 | sed 's:/$::g')
PR_2=$(echo $5 | sed 's:/$::g')

PEAK_SUFFIX=$(echo "${MERGED_PEAKS}" | perl -F'\.' -ane 'print $F[-1];')
OFPREFIX=$(basename "${MERGED_PEAKS}" | sed s/\."${PEAK_SUFFIX}"//g)

# Check for and make the necessary directories
IMD_FILES_DIR="consensusPeakCalling_intermediateFiles"
[[ ! -d "${IMD_FILES_DIR}" ]] && mkdir "${IMD_FILES_DIR}"

# =============================
# Keep any peaks that were called in all bioReps
# =============================
IN_ALL_BIO_REPS="${OFPREFIX}_inAllBioReps.${PEAK_SUFFIX}"
coverageBed -a ${BIO_REP_A} -b ${MERGED_PEAKS} | awk -v min=${OVERLAP_PORTION} '{if ($NF > min) print $0}' | coverageBed -a ${BIO_REP_B} -b - | awk -v min=${OVERLAP_PORTION} '{if ($NF > min) print $0}' | cut -f 1-3 > ${IN_ALL_BIO_REPS}

# =============================
# Remove peaks that don't intersect both the psuedoreplicates
# =============================
IN_BOTH_PRS="${OFPREFIX}_inBothPRs.${PEAK_SUFFIX}"
coverageBed -a ${PR_1} -b ${MERGED_PEAKS} | awk -v min=${OVERLAP_PORTION} '{if ($NF > min) print $0}' | coverageBed -a ${PR_2} -b - | awk -v min=${OVERLAP_PORTION} '{if ($NF > min) print $0}' |cut -f 1-3 > ${IN_BOTH_PRS}


# =============================
# Remove peaks from the psuedorep supported peaks that don't intersect at least one of the biological replicates
# And generate a txt file with a summary the intersections
# =============================
IN_A="${OFPREFIX}_inA.${PEAK_SUFFIX}"
IN_B="${OFPREFIX}_inB.${PEAK_SUFFIX}"
CONSENSUS_PEAKS="${OFPREFIX}_consensusPeaks.${PEAK_SUFFIX}"
PEAKS_LOST_TO_BLACKLIST="${OFPREFIX}_consensusPeaks_intersectingBlackList.${PEAK_SUFFIX}"
SUMMARY_FILE="${OFPREFIX}_intersectionSummary.txt"

# Get all peaks that are in both PRs and any 2 bio reps.  There should be a lot of duplication. I also only take the position information to avoid odd numbers of rows
coverageBed -a ${BIO_REP_A} -b ${IN_BOTH_PRS} | awk -v min=${OVERLAP_PORTION} '{if ($NF > min) print $0}' | cut -f 1-3 > ${IN_A}
coverageBed -a ${BIO_REP_B} -b ${IN_BOTH_PRS} |awk -v min=${OVERLAP_PORTION} '{if ($NF > min) print $0}' | cut -f 1-3 > ${IN_B}

# Get the original merged peaks (there summit info, etc) that intersect any of those listed in IN_2 (bedtools deals with the duplication)
cat ${IN_ALL_BIO_REPS} ${IN_A} ${IN_B} | bedtools sort -i - | intersectBed -wa -u -f 1 -a ${MERGED_PEAKS} -b - | intersectBed -wa -v -a - -b /Volumes/extra/Genomic_files/worm/ce10/ce10-blacklist.bed > ${CONSENSUS_PEAKS}
cat ${IN_ALL_BIO_REPS} ${IN_A} ${IN_B} | bedtools sort -i - | intersectBed -wa -u -f 1 -a ${MERGED_PEAKS} -b - | intersectBed -wa -u -a - -b /Volumes/extra/Genomic_files/worm/ce10/ce10-blacklist.bed > ${PEAKS_LOST_TO_BLACKLIST}

# write out details to the summary file
echo "Peaks in merged: " > ${SUMMARY_FILE}
cat ${MERGED_PEAKS} | wc -l >> ${SUMMARY_FILE}
echo "Peaks in both psuedoreps: " >> ${SUMMARY_FILE}
cat ${IN_BOTH_PRS} | wc -l >> ${SUMMARY_FILE}
echo "Peaks in both psuedoReps and bioRep A: " >> ${SUMMARY_FILE}
cat ${IN_A} | wc -l >> ${SUMMARY_FILE}
echo "Peaks in both psuedoReps and bioRep B: " >> ${SUMMARY_FILE}
cat ${IN_B} | wc -l >> ${SUMMARY_FILE}
echo "Peaks in all biological replicates: " >> ${SUMMARY_FILE}
cat ${IN_ALL_BIO_REPS} | wc -l >> ${SUMMARY_FILE}
echo "Consensus peaks lost to the blacklist: " >> ${SUMMARY_FILE}
cat ${PEAKS_LOST_TO_BLACKLIST} | wc -l >> ${SUMMARY_FILE}
echo "Peaks in consensus: " >> ${SUMMARY_FILE}
cat ${CONSENSUS_PEAKS} | wc -l >> ${SUMMARY_FILE}

mv ${PEAKS_LOST_TO_BLACKLIST} ${IN_A} ${IN_B} ${IN_BOTH_PRS} ${IN_ALL_BIO_REPS} ${IMD_FILES_DIR}
