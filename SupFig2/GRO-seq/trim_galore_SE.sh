#!/bin/bash
echo "Using trim_galore to filter paired FASTQ files" 1>&2

if [[ "$#" -lt 1 ]]
then
    echo "$(basename $0) [FASTQDir]"  1>&2
    echo "[FASTQDir]: directory containing FASTQ files to be filtered" 1>&2
    exit 1
fi

FASTQ_DIR=$(echo $1 | sed 's:/$::g')
cd ${FASTQ_DIR}

# make output directory if it doesnt exist
[[ ! -d "${FASTQ_DIR}/fastqc" ]] && mkdir "${FASTQ_DIR}/fastqc"


for f in $(find "$FASTQ_DIR" -name '*.fastq')
do
	echo "Processing: $f..."
	trim_galore --phred33 --quality 15 --stringency 15 --fastqc_args "--outdir ./fastqc"  $filePath/$f

done
