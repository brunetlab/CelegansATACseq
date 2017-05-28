#!/bin/bash

#### This will use tophat2 to map reads to C. elegans ce10.
# The gff3 used here is from gerstein et al, 2014
if [[ "$#" -lt 2 ]]
then
	echo "USEAGE: SE or PE read mapping:" 1>&2
	echo "This program accepts trimGalore'd fastq files and aligns the reads to ce10Mito with bowtie2 via tophat2" 1>&2
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

# 1) Map with tophat2

for FASTQ_FILE_1 in $(find "$FASTQ_DIR" -name '*_1_val_1.fq')
do
	OFPREFIX=$(basename "${FASTQ_FILE_1}" | sed 's/_1_val_1\.fq//g')
	FASTQ_FILE_2="${FASTQ_DIR}/${OFPREFIX}_2_val_2.fq"
	echo "Aligning: $OFPREFIX as Paired End"
	tophat2 -o "${BAM_DIR}/${OFPREFIX}" -p 3 --solexa-quals -g 10 --no-coverage-search --no-novel-indels --segment-length 17 --no-novel-juncs --b2-very-sensitive --transcriptome-index /Users/acd13/Softwares/tophat-2.0.9/transcriptomeIndexes/ce10_UcscRefSeqGenes_26Oct2014/UcscRefSeqGenes_26Oct2014 /Users/acd13/Softwares/bowtie2-2.1.0/indexes/ce10Mito ${FASTQ_FILE_1} ${FASTQ_FILE_2}
	
done