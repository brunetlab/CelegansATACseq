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

# Unfortunately I need to create a different set of binary files for each cell mark combination and then combine them afterwards
while read line ;do
	SHIFT=$(echo $line | perl -lane 'print $F[-1];')
	CELL=$(echo $line | perl -lane 'print $F[0];')
	MARK=$(echo $line | perl -lane 'print $F[1];')
	
	# Everything, but the first column
	echo $line | awk '{$NF=""; print $0}' | sed 's/ /	/g' > temp.txt 
	SAMPLE_ODIR="${ODIR}/${CELL}/${MARK}"  
	bash /Users/acd13/Desktop/ChromHMM/scriptsFromAnshul/myVersions/chmm_binarize.sh temp.txt ${SHIFT} ${SAMPLE_ODIR} /Volumes/extra/Genomic_files/worm/ce10/ce10.chrom.sizes 8 ${IDIR} ${POIS} ${IDIR}
done < ${FILE}
