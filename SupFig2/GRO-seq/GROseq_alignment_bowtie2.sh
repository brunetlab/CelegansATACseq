#!/bin/bash

# this takes the post-TrimGalore! reads then aligns the reads to ce10

if [[ "$#" -lt 2 ]]
then
	echo "USEAGE: SE read mapping:" 1>&2
	echo "This program accepts trimGalore'd fastq files and aligns the reads to ce10 with bowtie2" 1>&2
    echo "$(basename $0) [FASTQDir] [BAMDir]"  1>&2
    echo "   [FASTQDir]: directory containing FASTQ files to be aligned; (with format: *_1.trimmed.fq)" 1>&2
    echo "   [BAMDir]: directory to put BAM files in" 1>&2
    exit 1
fi

FASTQ_DIR=$(echo $1 | sed 's:/$::g')
BAM_DIR=$(echo $2 | sed 's:/$::g')

# make bam directory if it doesnt exist
[[ ! -d "${BAM_DIR}/" ]] && mkdir "${BAM_DIR}/"
[[ ! -d "${BAM_DIR}/metrics" ]] && mkdir "${BAM_DIR}/metrics"

# 1) Map with bowite2

cp /Users/acd13/Softwares/bowtie2-2.1.0/indexes/ce10mito*.bt2 ./
for f in $(find "$FASTQ_DIR" -name '*_1_trimmed.fq')
do
	echo "Aligning: $f...\n"
	FASTQ_FILE_1=$f
	OFPREFIX=$(basename "${f}" | sed 's/_1\.trim\.fastq//g')
	
	RAW_SAM_FILE="${OFPREFIX}.raw.sam"
	ALIGNMENT_METRICS_FILE="${OFPREFIX}.alignMetrics.txt"
	
	bowtie2 -p 4 --very-sensitive --met-file ${ALIGNMENT_METRICS_FILE} -q -x ce10mito -U ${FASTQ_FILE_1} -S ${RAW_SAM_FILE}
	
	# ==============================================================                                                                                                                                                                                                              
	# Sort by position                                                                                                                                                                                                               
	# ============================================================== 
	RAW_BAM_PREFIX="${OFPREFIX}.raw.srt"
	RAW_BAM_FILE="${RAW_BAM_PREFIX}.bam" # To be stored                                                                                                                                                                                                                           
	RAW_BAM_FILE_MAPSTATS="${RAW_BAM_PREFIX}.flagstat.qc" # QC File

	samtools view -Sb ${RAW_SAM_FILE} | samtools sort - ${RAW_BAM_PREFIX}
	
	rm ${RAW_SAM_FILE}

	samtools flagstat ${RAW_BAM_FILE} > ${RAW_BAM_FILE_MAPSTATS}
	
	mv ${RAW_BAM_FILE} ${BAM_DIR}
	mv ${RAW_BAM_FILE_MAPSTATS} ${ALIGNMENT_METRICS_FILE} "${BAM_DIR}/metrics"

done

rm ce10mito*.bt2
echo "Finished mapping"
