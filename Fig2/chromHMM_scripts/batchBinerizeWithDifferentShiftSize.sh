#!/bin/bash

if [[ "$#" -lt 3 ]]
then
    echo "USAGE: $(basename $0) [METADATA] [poissonThreshold] [ODIR]" 1>&2
    echo 'Will take a chromHMM metadata file with an extra column: the shift size, and run the binerize command using 8GB of memory' 1>&2
    echo "$(basename $0) [METADATA] [ODIR]"  1>&2
    echo "   [METADATA]: Txt file containing the metaData for binerizing (see chromHMM manual)," 1>&2
    echo "   plus an additional last column: the shift size (frag length/2)." 1>&2
    echo '   [poissonThreshold]:  specifies the tail probability of the poisson distribution that the binarization threshold should correspond to. The default value of this parameter is 0.0001' >&2
    echo "   [ODIR]: Output directory" 1>&2
    exit 1
fi

FILE=$(echo $1 | sed 's:/$::g')
POIS=$(echo $2 | sed 's:/$::g')
ODIR=$(echo $3 | sed 's:/$::g')

[[ ! -d ${ODIR} ]] && mkdir -p ${ODIR}



IDIR=$(dirname ${FILE})
cd ${IDIR}

while read line ;do
	SHIFT=$(echo $line | perl -lane 'print $F[-1];')
	echo $line | awk '{$NF=""; print $0}' | sed 's/ /	/g' > temp.txt 
	TEMP_ODIR="${ODIR}/temp"  
	bash /Users/acd13/Desktop/ChromHMM/scriptsFromAnshul/myVersions/chmm_binarize.sh temp.txt ${SHIFT} ${TEMP_ODIR} /Volumes/extra/Genomic_files/worm/ce10/ce10.chrom.sizes 8 ${IDIR} ${POIS} ${IDIR}
   
	for BIN_FILE in $(find "$TEMP_ODIR" -name '*binary.txt')
	do
		FILENAME=$(basename "${BIN_FILE}")
		EXISTING_FILE="${ODIR}/${FILENAME}"
		if [[ -f ${EXISTING_FILE} ]] 
		then 
			paste ${EXISTING_FILE} ${BIN_FILE} > "${ODIR}/tmp.txt"
	    	mv "${ODIR}/tmp.txt" ${EXISTING_FILE}
	    else
	    	mv ${BIN_FILE} ${EXISTING_FILE}
	    fi
	    
	done

done < ${FILE}

rm -r temp.txt "${ODIR}/temp"

# the problem with pasting is that the first line gets all messed up, it should only contain cell type and chromosome once, so now we fix it
for BIN_FILE in $(find "$ODIR" -name '*binary.txt')
	do
		NEW_HEADER=$(cut -f 1-2 ${BIN_FILE} | head -1) # get the right header
		echo $NEW_HEADER > tempBin.txt
		tail -n +2 ${BIN_FILE} >> tempBin.txt
		
		mv tempBin.txt ${BIN_FILE}
	    
	done