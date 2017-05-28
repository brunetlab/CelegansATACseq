#!/bin/bash

#### This takes a directory of PE fastq ATAC reads from the sequencer,
# trims adapters, aligns, filters for quality alignment,
# generates insert sizes (a good metric for the quality of the library prep),
# adjusts for the Tn5 binding footprint,
# and creates bed and bam files of just the insert site (i.e. the 5' end of each read).
# Finally a bed file of inert sites shifted 75bp for use with MACS2 to call peaks is created
# Peaks are not called however.  To do that you should take the output from this and submit replicates to ___.sh

if [[ "$#" -lt 2 ]]
then
	echo "USEAGE: ATAC-seq read trimming, mapping, filtering, adjusting, trimming to insert sites and creating files for subsequent peak calling with MACS2:" 1>&2
	echo "This program accepts the original PE fastq files and takes care of evrything else up until peak calling." 1>&2
    echo "$(basename $0) [FASTQDir] [BAMDir]"  1>&2
    echo "   [FASTQDir]: directory containing FASTQ files to be trimmed and aligned; (with format: *_1_pf.fastq and *_2_pf.fastq)" 1>&2
    echo "   [BAMDir]: directory to put BAM files in" 1>&2
    exit 1
fi

# =============================
# These need to be customized for each users environment
# ==============================
BOWTIE2_INDEXES="/Users/acd13/Softwares/bowtie2-2.1.0/indexes/ce10mito"
MAPQ_THRESH=30 # This is the minimum MAPQ score to be retained, I wouldn't recommend changing this unless you have a definite reason.
PICARD_PATH="/Users/acd13/Softwares/picard-tools-1.74/picard-tools-1.74"
GENOME_SIZE="/Volumes/extra/Genomic_files/worm/ce10/ce10.chrom.sizes"
GENOME_NAME="ce10"
 
# =============================
# Read in directories and create what's needed
# ==============================
FASTQ_DIR=$(echo $1 | sed 's:/$::g')
BAM_DIR=$(echo $2 | sed 's:/$::g')
INSERT_DIR="${BAM_DIR}/SingleBPInserts"
# make bam directory if it doesnt exist
[[ ! -d "${FASTQ_DIR}/Raw" ]] && mkdir "${FASTQ_DIR}/Raw"
[[ ! -d "${BAM_DIR}/" ]] && mkdir "${BAM_DIR}/"
[[ ! -d "${BAM_DIR}/metrics" ]] && mkdir "${BAM_DIR}/metrics"
[[ ! -d "${BAM_DIR}/InsertSizeMetrics" ]] && mkdir "${BAM_DIR}/InsertSizeMetrics"
[[ ! -d "${INSERT_DIR}/" ]] && mkdir "${INSERT_DIR}/"
[[ ! -d "${INSERT_DIR}/HomerTags" ]] && mkdir "${INSERT_DIR}/HomerTags"

# Set up this too
ATAC_SCRIPTS_PATH="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" # that is, the directory this script is saved in

#================================
# Define functions
#================================

function getInsertSites {	
	# Files and directories to write to
	INSERT_SITES_BED_FILE="${2}.shifted.insertSites.bed"
	INSERT_SITES_BAM_PREFIX="${2}.shifted.insertSites.sort"
	INSERT_SITES_BAM_FILE="${2}.shifted.insertSites.sort.bam"
	HOMER_TAG_DIR="${INSERT_DIR}/HomerTags/${2}"
	
	# Convert the bam to bed, but a a custom function that lists the insert size as the score
	# Then adjust 4 or 5bp for the Tn5 binding site, and finally trim to the 5prime end of the read.
	bamToBedInsertLengthAsScore "${1}" | adjustBedTn5 | FivePrimOfBed6Cols > "${INSERT_SITES_BED_FILE}"
	# Convert back to a bam, but now just the 1bp
	bedToBam -i "${INSERT_SITES_BED_FILE}" -g "${GENOME_SIZE}" > temp.bam
	samtools sort temp.bam "${INSERT_SITES_BAM_PREFIX}"
	rm temp.bam
	samtools index "${INSERT_SITES_BAM_FILE}"
	# This can be commented out if you don't want a Homer tag directory.
	makeTagDirectory "${HOMER_TAG_DIR}" -genome "${GENOME_NAME}" -keepAll -checkGC -fragLength 1 "${INSERT_SITES_BAM_FILE}"
	
	# And then a little clean up to save space
	gzip "${INSERT_SITES_BED_FILE}"
}

function bamToBedInsertLengthAsScore {
	bamToBed -i "${1}"> temp.bed 
	samtools view  "${1}" | perl -lane 'print $F[0]."\t".$F[8];' > insert.temp
	
	paste temp.bed insert.temp | \
	perl -lane '$name=(split/\//,$F[3])[0]; if ($name =~ /$F[6]/){ print $F[0]."\t".$F[1]."\t".$F[2]."\t".$F[3]."\t".$F[7]."\t".$F[5]; }else{ print STDERR "NAMES DONT MATCH for $name $F[6]!!!!";}'
	rm temp.bed insert.temp
}

function adjustBedTn5 {
	awk -F $'\t' 'BEGIN {OFS = FS}{if ($6 == "+") { $2 = $2 + 4 } else if ($6 == "-") {$3 = $3 - 5} print $0}' "${1}"
}

function fivePrimeOfBed6Cols {
	cat "${1}" | awk 'BEGIN{OFS="\t"}{if($6 == "-") $2=$3-1; print $1, $2, $2+1, $4, $5, $6}' | sortBed -i -
}

# find bad CIGAR read names, and remove them, then sort the resulting bam
function removeReadsWithBadCigars {

	cat "${1}" | awk 'BEGIN {FS="\t" ; OFS="\t"} ! /^@/ && $6!="*" { cigar=$6; gsub("[0-9]+D","",cigar); n = split(cigar,vals,"[A-Z]"); s = 0; for (i=1;i<=n;i++) s=s+vals[i]; seqlen=length($10) ; if (s!=seqlen) print $1 ; }' | sort | uniq > badCigars.temp

	# Remove bad CIGAR read pairs
	if [[ $(cat badCigars.temp | wc -l) -gt 0 ]]
	then
		cat "${1}" | grep -v -F -f badCigars.temp | samtools view -Sb - | samtools sort - "${2}"
	else
		samtools view -Sb "${1}" | samtools sort - "${2}"
	fi

	rm badCigars.temp "${1}"
	
}

# =============================
# Trim adaptor sequences
# ==============================
TRIM_ADAPTORS="${ATAC_SCRIPTS_PATH}/pyadapter_trim.py"
cd "${FASTQ_DIR}"
for FIRST_FQ in $(find "${FASTQ_DIR}" -name '*_1_pf.fastq')
do
	# Get file names
	FILE_PREFIX=$(basename "${FIRST_FQ}" | sed 's/_1_pf\.fastq//g')
	SECOND_FQ="${FILE_PREFIX}_2_pf.fastq"
	
	# Trim the adaptors
	echo "Trimming adaptors from: ${FILE_PREFIX}"
	python "${TRIM_ADAPTORS}" -a "${FIRST_FQ}" -b "${SECOND_FQ}"
	
	# zip the raw files and move them to a separate directory
	SECOND_GZ="${FILE_PREFIX}_2_pf.fastq.gz"
	FIRST_GZ="${FILE_PREFIX}_1_pf.fastq.gz"
	gzip "${SECOND_FQ}"
	gzip "${FIRST_FQ}"
	mv "${FIRST_GZ}" "${FASTQ_DIR}/Raw/"
	mv "${SECOND_GZ}" "${FASTQ_DIR}/Raw/"
done
echo "Finished trimming"

# =============================
# Map with bowite2
# ==============================
for FASTQ_FILE_1 in $(find "${FASTQ_DIR}" -name '*_1_pf.trim.fastq')
do
	FASTQ_FILE_2=$(basename "${FASTQ_FILE_1}" | sed 's/_1_pf\.trim\.fastq/_2_pf\.trim\.fastq/g')
	OFPREFIX=$(basename "${FASTQ_FILE_1}" | sed 's/_1_pf\.trim\.fastq//g')
	
	echo "Aligning: ${OFPREFIX}"
	
	RAW_SAM_FILE="${OFPREFIX}.sam"
	ALIGNMENT_METRICS_FILE="${OFPREFIX}.alignMetrics.txt"
	
	bowtie2 -p 4 --end-to-end --no-mixed -X 2000 --met 1 --met-file "${ALIGNMENT_METRICS_FILE}" -q -x "${BOWTIE2_INDEXES}" -1 "${FASTQ_FILE_1}" -2 "${FASTQ_FILE_2}" -S "${RAW_SAM_FILE}"
	
	#=============================================================                                                                                                                 
	# Remove read pairs with bad CIGAR strings and sort by position
	# ============================================================== 
	RAW_BAM_PREFIX="${OFPREFIX}"
	RAW_BAM_FILE="${RAW_BAM_PREFIX}/Users/acd13/Desktop/ATAC/bowtie2Bams/N2_dev/thirdRep"
	BADCIGAR_FILE="${RAW_BAM_PREFIX}.badReads.tmp"
	RAW_BAM_FILE_MAPSTATS="${RAW_BAM_PREFIX}.flagstat.qc" # QC File
	
	removeReadsWithBadCigars "${RAW_SAM_FILE}" "${RAW_BAM_PREFIX}"
	samtools flagstat "${RAW_BAM_FILE}" > "${RAW_BAM_FILE_MAPSTATS}"
	
	#==============================================================                                                                                                                
	# Clean up FASTQ dir
	# ==============================================================
	
	mv "${RAW_BAM_FILE}" "${BAM_DIR}"
	mv "${RAW_BAM_FILE_MAPSTATS}" "${ALIGNMENT_METRICS_FILE}" "${BAM_DIR}/metrics"
	gzip "${FASTQ_FILE_1}"
	gzip "${FASTQ_FILE_2}"

done
	

#==============================================================                                                                                                                
# Fully process the bams
# ==============================================================

cd "${BAM_DIR}"

for RAW_BAM_FILE in $(find "${BAM_DIR}" -name '*.bam')
do

	OFPREFIX=$(basename "${RAW_BAM_FILE}" | sed 's/\.bam//g')    

    # =============================
    # Remove  unmapped, mate unmapped
    # not primary alignment, reads failing platform
    # Remove low MAPQ reads
    # Obtain name sorted BAM file
    # ==================
    FILT_BAM_PREFIX="${OFPREFIX}.filt" 
    FILT_BAM_FILE="${FILT_BAM_PREFIX}.bam"

    samtools view -F 1804 -q "${MAPQ_THRESH}" -b "${RAW_BAM_FILE}" | samtools sort - "${FILT_BAM_PREFIX}"

    # ========================
    # Mark duplicates
    # ======================
    TMP_FILT_BAM_FILE="${FILT_BAM_PREFIX}.dupmark.bam"
    DUP_FILE_QC="${FILT_BAM_PREFIX}.dup.qc"

    java -Xmx4G -jar "$PICARD_PATH/MarkDuplicates.jar" INPUT="${FILT_BAM_FILE}" OUTPUT="${TMP_FILT_BAM_FILE}" METRICS_FILE="${DUP_FILE_QC}" VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false

    mv "${TMP_FILT_BAM_FILE}" "${FILT_BAM_FILE}"
    
    # ============================
    # Remove duplicates
    # Index final position sorted BAM
    # ============================
    FINAL_BAM_PREFIX="${OFPREFIX}.filt.nodup"
    FINAL_BAM_FILE="${FINAL_BAM_PREFIX}.bam"
    FINAL_BAM_INDEX_FILE="${FINAL_BAM_PREFIX}.bai"
    FINAL_BAM_FILE_MAPSTATS="${FINAL_BAM_PREFIX}.flagstat.qc"

    samtools view -F 1804 -b "${FILT_BAM_FILE}" > "${FINAL_BAM_FILE}"

    # Index Final BAM file
    samtools index "${FINAL_BAM_FILE}" "${FINAL_BAM_INDEX_FILE}"

	# generate the final stats
    samtools flagstat "${FINAL_BAM_FILE}" > "${FINAL_BAM_FILE_MAPSTATS}"
    
    # ============================
    # Get Insert Size metrics
    # ============================
    INSERT_TXT_FILE="${FINAL_BAM_PREFIX}.insertSizes.txt"
	INSERT_HISTO_FILE="${FINAL_BAM_PREFIX}.insertSizes.pdf"
	java -Xms1g -Xmx4g -jar "$PICARD_PATH/CollectInsertSizeMetrics.jar" METRIC_ACCUMULATION_LEVEL=ALL_READS OUTPUT="${INSERT_TXT_FILE}" HISTOGRAM_FILE="${INSERT_HISTO_FILE}" INPUT="${FINAL_BAM_FILE}"
    
    # ============================
    # Get insert sites - as a gzipped bed and a bam file as well as Homer tag directory
    # ============================
    INSERTS_PREFIX="${INSERT_DIR}/${FINAL_BAM_PREFIX}.adjusted.insertSites"
    INSERTS_BED_FILE="${INSERTS_PREFIX}.bed.gz"
    
	getInsertSites "${FINAL_BAM_FILE}" "${INSERTS_PREFIX}"

    # ============================
    # Get ready for calling peaks
    # ============================
    BED_FILE_FOR_MACS="${INSERTS_PREFIX}.75bpShiftedForMACS.bed"    
    
    zcat "${INSERTS_BED_FILE}" | slopBed -i - -g "${GENOME_SIZE}" -l 75 -r -75 -s > "${BED_FILE_FOR_MACS}"
    
    # ============================
    # Clean Up
    # ============================
    
    # move files to correct directories
    mv "${INSERT_TXT_FILE}" "${INSERT_HISTO_FILE}" "${BAM_DIR}/InsertSizeMetrics"
    mv "${DUP_FILE_QC}" "${FINAL_BAM_FILE_MAPSTATS}" "${BAM_DIR}/metrics"

done
