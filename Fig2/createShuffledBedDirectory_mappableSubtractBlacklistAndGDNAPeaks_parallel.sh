#!/Users/acd13/bin/parallel --shebang-wrap bash

# #!/bin/bash

#### Created 07.14.2014 - ACD
# this program will take a bedfile of peaks and create a folder of a specified number of gzipped shuffled beds to be used for enrichment calculations

if [[ "$#" -lt 2 ]]
then
	echo "USEAGE: Generate a folder of gzipped shuffled beds to be used for example for enrichment calculations:" 1>&2
	echo "This program accepts a bed file, the number of shuffled beds to create, and the output directory" 1>&2
    echo "$(basename $0) [BED_FILE]"  1>&2
    echo "$(basename $1) [OUTPUT_DIR]"  1>&2
    echo "   [BED_FILE]: bed file to shuffle" 1>&2
    echo "   [OUTPUT_DIR]: where to store the output" 1>&2
    exit 1
fi


# Get the arguments
BED_FILE=$(echo $1 | sed 's:/$::g')
OUTPUT_DIR=$(echo $2 | sed 's:/$::g')

PREFIX=$(basename "${BED_FILE}" | sed 's/\.bed//g')
TRIMMED_BED="${PREFIX}_trimmed.bed"

# Trim the bed file to the minimal info (i.e. just position)
cut -f 1-3 ${BED_FILE} > ${TRIMMED_BED}
	
function runShufBed {
	# Set the variables
	BED_FILE=$1
	OUTPUT_DIR=$2
	i=$3
	PREFIX=$(basename "${BED_FILE}" | sed 's/_trimmed\.bed//g')
	OUT_FILE="${OUTPUT_DIR}/${PREFIX}_shuffled.${i}.bed.gz"
	MAPPABLE_GENOME="/Volumes/extra/Genomic_files/worm/ce10/ce10.K100.mappable.subtract_blackList_N2GdnaMacs2.1q5e-2peaks0.5Overlap.bed"
	GENOME="/Volumes/extra/Genomic_files/worm/ce10/ce10.chrom.sizes"
		
	shuffleBed -incl ${MAPPABLE_GENOME} -i ${BED_FILE} -g ${GENOME} -noOverlapping -seed ${i} -maxTries 10000 | gzip -c > ${OUT_FILE}
	
}

#for i in {1..10000}
#do
#	OUT_FILE="${OUTPUT_DIR}/${OUT_FILE_PREFIX}.${i}.bed.gz"
#    shuffleBed -incl ${MAPPABLE_GENOME} -i ${TRIMMED_BED} -g ${GENOME} -noOverlapping -seed ${i} | gzip -c > ${OUT_FILE}
#done
  
export -f runShufBed
parallel 'runShufBed {}' ::: ${TRIMMED_BED} ::: ${OUTPUT_DIR} ::: {1..10000}
    
rm ${TRIMMED_BED}