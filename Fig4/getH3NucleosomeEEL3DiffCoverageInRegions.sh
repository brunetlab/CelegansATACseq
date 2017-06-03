#nucStage=$1 # must be EE or L3
regionOfInterest=$1 # must be a bed file

if [[ $regionOfInterest =~ \.gz$ ]]
then
	temp="${regionOfInterest#*.}"
	regionOfInterestSuffix=$(echo ${temp} | sed 's/\.gz//g')
	outFilePrefix=$(basename $regionOfInterest | sed 's/\.${regionOfInterestSuffix}\.gz//g')
else
	regionOfInterestSuffix="${regionOfInterest##*.}"
	outFilePrefix=$(basename $regionOfInterest | sed 's/\.${regionOfInterestSuffix}//g')
fi

echo "File ending is $regionOfInterestSuffix"

outputDir=$2 # need not exist

nucSignalDir='/Users/acd13/Desktop/ATAC/Analysis/nucleosomeCalls/H3_EE_v_L3_fromDanpos/result/diff/'
cd ${nucSignalDir}
nucSignalBed=EE_H3-L3_H3.pois_diff.wig.srt.bg

mkdir -p ${outputDir}

bedtools map -a ${regionOfInterest} -b ${nucSignalBed} -c 4 -o median > ${outputDir}/${outFilePrefix}_EEL3Diff_H3InferredNucleosomeMedianSignal.${regionOfInterestSuffix}
