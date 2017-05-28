#!/bin/bash

# This program takes raw bam files (e.g. from bowtie2) that have been sorted, and removes low quality reads, any that don't map, and then splits them strand specifically
# I'm also checking the duplicates after that, but we'll see if I remove those

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

    # =============================
    # Remove  unmapped, mate unmapped
    # not primary alignment, reads failing platform
    # Remove low MAPQ reads
    # Obtain name sorted BAM file
    # ==================
    ALL_FILT_BAM_PREFIX="${OFPREFIX}.filt.srt" # <name>.filt.srt.posStrand
    POS_FILT_BAM_PREFIX="${OFPREFIX}.filt.srt.posStrand" # <name>.filt.srt.posStrand
    NEG_FILT_BAM_PREFIX="${OFPREFIX}.filt.srt.negStrand" # <name>.filt.srt.negStrand
    ALL_FILT_BAM_FILE="${ALL_FILT_BAM_PREFIX}.bam"  # <name>.filt.srt.bam
    POS_FILT_BAM_FILE="${POS_FILT_BAM_PREFIX}.bam"  # <name>.filt.srt.posStrand.bam
    NEG_FILT_BAM_FILE="${NEG_FILT_BAM_PREFIX}.bam"  # <name>.filt.srt.negStrand.bam
    MAPQ_THRESH=30
    
    samtools view -q ${MAPQ_THRESH} -F 16 -b ${RAW_BAM_FILE}| samtools sort - ${POS_FILT_BAM_PREFIX}
    samtools view -q ${MAPQ_THRESH} -f 16 -b ${RAW_BAM_FILE}| samtools sort - ${NEG_FILT_BAM_PREFIX}
    samtools view -q ${MAPQ_THRESH} -b ${RAW_BAM_FILE}| samtools sort - ${ALL_FILT_BAM_PREFIX}

    
    
    # ========================
    # Mark duplicates - Positive strand, then negative strand
    # ========================
    
    TMP_FILT_BAM_FILE="${POS_FILT_BAM_PREFIX}.dupmark.bam" # <name>.filt.srt.dupmark.bam
    MARKDUP="/Users/acd13/Softwares/picard-tools-1.74/picard-tools-1.74/MarkDuplicates.jar"
    DUP_FILE_QC="${POS_FILT_BAM_PREFIX}.dup.qc" # QC file <name>.filt.srt.sup.qc

    java -Xmx4G -jar ${MARKDUP} INPUT=${POS_FILT_BAM_FILE} OUTPUT=${TMP_FILT_BAM_FILE} METRICS_FILE=${DUP_FILE_QC} VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false

    mv ${TMP_FILT_BAM_FILE} ${POS_FILT_BAM_FILE}

	# ============================
	# Remove duplicates
	# Index final position sorted BAM
	# ============================
	
	FINAL_BAM_PREFIX="${POS_FILT_BAM_PREFIX}.nodup" 
	FINAL_BAM_FILE="${FINAL_BAM_PREFIX}.bam" 
	FINAL_BAM_INDEX_FILE="${FINAL_BAM_PREFIX}.bai" 
	FINAL_BAM_FILE_MAPSTATS="${FINAL_BAM_PREFIX}.flagstat.qc" 

	samtools view -F 1804 -b ${POS_FILT_BAM_FILE} > ${FINAL_BAM_FILE}

	# Index Final BAM file
	samtools index ${FINAL_BAM_FILE} ${FINAL_BAM_INDEX_FILE}

	samtools flagstat ${FINAL_BAM_FILE} > ${FINAL_BAM_FILE_MAPSTATS}
	
	
	 # Now negative strand
	TMP_FILT_BAM_FILE="${NEG_FILT_BAM_PREFIX}.dupmark.bam" # <name>.filt.srt.dupmark.bam
    DUP_FILE_QC="${NEG_FILT_BAM_PREFIX}.dup.qc" # QC file <name>.filt.srt.sup.qc

    java -Xmx4G -jar ${MARKDUP} INPUT=${NEG_FILT_BAM_FILE} OUTPUT=${TMP_FILT_BAM_FILE} METRICS_FILE=${DUP_FILE_QC} VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false

    mv ${TMP_FILT_BAM_FILE} ${NEG_FILT_BAM_FILE}

	# ============================
	# Remove duplicates
	# Index final position sorted BAM
	# ============================
	
	FINAL_BAM_PREFIX="${NEG_FILT_BAM_PREFIX}.nodup" 
	FINAL_BAM_FILE="${FINAL_BAM_PREFIX}.bam" 
	FINAL_BAM_INDEX_FILE="${FINAL_BAM_PREFIX}.bai" 
	FINAL_BAM_FILE_MAPSTATS="${FINAL_BAM_PREFIX}.flagstat.qc" 

	samtools view -F 1804 -b ${NEG_FILT_BAM_FILE} > ${FINAL_BAM_FILE}

	# Index Final BAM file
	samtools index ${FINAL_BAM_FILE} ${FINAL_BAM_INDEX_FILE}

	samtools flagstat ${FINAL_BAM_FILE} > ${FINAL_BAM_FILE_MAPSTATS}
	
	

done

