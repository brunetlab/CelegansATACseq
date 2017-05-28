#!/bin/bash

# This program takes bam files and removes duplicates

if [[ "$#" -lt 1 ]]
then
    echo "$(basename $0) [BAM_DIR]"  1>&2
    echo "   [BAM_DIR]: directory containing BAM files to be processed" 1>&2
    exit 1
fi

BAM_DIR=$(echo $1 | sed 's:/$::g')

# 1) Filter bams to only mapped reads, remove duplicates and calculate some QCs
cd ${BAM_DIR}

for f in $(find "$BAM_DIR" -name '*.bam')
do
    OFPREFIX=$(basename "${f}" | sed 's/\.bam//g')
	RAW_BAM_FILE="${f}"
    
    # ========================
    # Mark duplicates 
    # ========================
    
    TMP_FILT_BAM_FILE="${OFPREFIX}.dupmark.bam" 
    MARKDUP="/Users/acd13/Softwares/picard-tools-1.74/picard-tools-1.74/MarkDuplicates.jar"
    DUP_FILE_QC="${OFPREFIX}.dup.qc" # QC file <name>.filt.srt.sup.qc

    java -Xmx4G -jar ${MARKDUP} INPUT=${RAW_BAM_FILE} OUTPUT=${TMP_FILT_BAM_FILE} METRICS_FILE=${DUP_FILE_QC} VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false

    mv ${TMP_FILT_BAM_FILE} ${RAW_BAM_FILE}

	# ============================
	# Remove duplicates
	# Index final position sorted BAM
	# ============================
	
	FINAL_BAM_PREFIX="${OFPREFIX}.nodup" 
	FINAL_BAM_FILE="${FINAL_BAM_PREFIX}.bam" 
	FINAL_BAM_INDEX_FILE="${FINAL_BAM_PREFIX}.bai" 
	FINAL_BAM_FILE_MAPSTATS="${FINAL_BAM_PREFIX}.flagstat.qc" 

	samtools view -F 1804 -b ${RAW_BAM_FILE} > ${FINAL_BAM_FILE}

	# Index Final BAM file
	samtools index ${FINAL_BAM_FILE} ${FINAL_BAM_INDEX_FILE}

	samtools flagstat ${FINAL_BAM_FILE} > ${FINAL_BAM_FILE_MAPSTATS}
	

done

