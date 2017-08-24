#!/bin/bash

OVERLAP_PORTION="0.5" # that is the portion of a merged peak that must be intersected by a peak in the comparison files (either the reps or the prs)
MACS_SUFFIX="MACS2.1_q5e-2" # This is what is added to the files created by MACS2.1
MACS2_QVAL="5e-2"
SCRIPTS_PATH="/Users/acd13/Dropbox/Scripts/ATAC_scripts"
TOTAL_SIZE='9e7' # The size of the useable genome for MACS2 this is for ce10
GDNA_OVERLAP_PORTION="0.2" # that is the portion of a consensus peak that must be intersected by a gDNA peak to remove that consensus peak
BLACKLIST_PEAKS="/Volumes/extra/Genomic_files/worm/ce10/ce10-blacklist.bed"

function callSubPeaks {
	BED_FILE=$(echo "${1}" | sed 's:/$::g')
	SUFFIX=$(echo "${2}" | sed 's:/$::g')
	PREFIX=$(basename "${BED_FILE}" | sed 's/\.bed//g')

	MACS_OUTPUT_PREFIX="${PREFIX}.${SUFFIX}"
    macs2 callpeak -t "${BED_FILE}" -f BED -n "${MACS_OUTPUT_PREFIX}" -g "${TOTAL_SIZE}" -q "${MACS2_QVAL}" --nomodel --extsize 150 -B --keep-dup all --call-summits
    
	# ============================
	# Process Peaks
	# ============================
	MACS_PEAKS="${MACS_OUTPUT_PREFIX}_peaks.narrowPeak"
	SUBPEAKS_SPLITTER="${SCRIPTS_PATH}/splitMACS2SubPeaks.pl" # this code simply takes the mid point between two subpeaks
	SPLIT_PEAKS="${MACS_OUTPUT_PREFIX}_splitSubsPeaks.narrowPeak"
	BLACKLIST_CLEARED_PEAKS="${MACS_OUTPUT_PREFIX}_splitSubsPeaks.blacklistCleared.narrowPeak"
	BLACKLIST_LOST_PEAKS="${MACS_OUTPUT_PREFIX}_splitSubsPeaks.blacklistLost.narrowPeak"

	perl "${SUBPEAKS_SPLITTER}" "${MACS_PEAKS}" > "${SPLIT_PEAKS}"
	
	intersectBed -v -a "${SPLIT_PEAKS}" -b "${BLACKLIST_PEAKS}" > "${BLACKLIST_CLEARED_PEAKS}"
	intersectBed -u -a "${SPLIT_PEAKS}" -b "${BLACKLIST_PEAKS}" > "${BLACKLIST_LOST_PEAKS}"
 
}

MERGED_PREFIX="allReps"
# I previously did this
MERGED_REPS="/Users/acd13/Desktop/ATAC/insertSites/N2dev/shiftedForMacs2/chromEndsRemoved/allReps.adjusted.insertSites.75bpShiftedForMACS.bed"
#cat *.filt.nodup.adjusted.insertSites.75bpShiftedForMACS.bed | sort -k 1,1 -k2,2n > "${MERGED_REPS}"

PR_1="${MERGED_PREFIX}_pr1.bed"
PR_2="${MERGED_PREFIX}_pr2.bed"
perl /Users/acd13/Dropbox/Scripts/MostUsed/random_split_bed.pl "${MERGED_REPS}" "${PR_1}" "${PR_2}"

MERGED_REPS_PEAKS="/Users/acd13/Desktop/ATAC/macs2Peaks/usingInserts_75bpShift/N2_dev/withAllReps/allReps/allReps.MACS2.1_q5e-2_peaks.splitSubsPeaks.blacklistCleared.narrowPeak"
PR1_PREFIX=$(basename "${PR_1}" | sed 's/\.bed//g')
PR1_PEAKS="${PR1_PREFIX}.${MACS_SUFFIX}_splitSubsPeaks.blacklistCleared.narrowPeak"
PR2_PREFIX=$(basename "${PR_1}" | sed 's/\.bed//g')
PR2_PEAKS="${PR2_PREFIX}.${MACS_SUFFIX}_splitSubsPeaks.blacklistCleared.narrowPeak"

for BED_FILE in "${PR_1}" "${PR_2}" #"${MERGED_REPS}" I had previously created this
do
	callSubPeaks ${BED_FILE} ${MACS_SUFFIX}
done

# ============================
# Now call the consensus peaks
# And remove any of those intersecting gDNA peaks
# ============================
# The final product is this:
CONSENSUS_PREFIX="${MERGED_PREFIX}.${MACS_SUFFIX}_splitSubsPeaks_blacklistCleared_consensusPeaks"
CONSENSUS_PEAKS="${CONSENSUS_PREFIX}.narrowPeak"
coverageBed -a "${PR1_PEAKS}" -b "${MERGED_REPS_PEAKS}" | awk -v min="${OVERLAP_PORTION}" '{if ($NF > min) print $0}' | \
coverageBed -a "${PR2_PEAKS}" -b - | awk -v min="${OVERLAP_PORTION}" '{if ($NF > min) print $0}' |\
cut -f 1-3 > "${CONSENSUS_PEAKS}"

# Now remove gDNA peaks
GDNA_PEAKS="/Users/acd13/Desktop/ATAC/macs2Peaks/usingInserts_75bpShift/N2_dev/withAllReps/gDNAMasked_MACS5e-2_0.5Overlap/N2_gDNA.filt.nodup.adjusted.macs2.1q5e-2_peaks.narrowPeak"
FINAL_PEAKS="${CONSENSUS_PREFIX}_gDNA${GDNA_OVERLAP_PORTION}PeaksRemoved.narrowPeak"
GDNA_LOST_PEAKS="${CONSENSUS_PREFIX}_LostToGDNA${GDNA_OVERLAP_PORTION}Peaks.narrowPeak"
intersectBed -wa -v -f "${GDNA_OVERLAP_PORTION}" -a "${CONSENSUS_PEAKS}" -b "${GDNA_PEAKS}" > "${FINAL_PEAKS}"
intersectBed -wa -u -f "${GDNA_OVERLAP_PORTION}" -a "${CONSENSUS_PEAKS}" -b "${GDNA_PEAKS}" > "${GDNA_LOST_PEAKS}"

