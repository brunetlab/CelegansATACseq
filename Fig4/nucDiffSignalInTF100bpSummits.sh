cd /Users/acd13/Desktop/ATAC/Analysis/nucleosomeCalls
mkdir -p nucleosomeSignalInTFPeaks
cd nucleosomeSignalInTFPeaks

tfPeakDir='/Users/acd13/Desktop/ATAC/Analysis/tfs/qualityTfs/peaksToUse/'
for peakStage in L3 EE YA
do
    mkdir -p ${peakStage}TF100bpSummits
   	for f in ${tfPeakDir}${peakStage}/*.narrowPeak.gz
   	do
   		tempFile=$f.100bpsummit.narrowPeak
 		# b/c these are narrowPeaks, they have the summits as the last column
       	zcat ${f} | perl -lane '$newStart = ($F[1] + $F[-1] - 50); print join("\t", ($F[0], $newStart, ($newStart + 100)));' | sortBed -i - > ${tempFile}
        bash ../getH3NucleosomeEEL3DiffCoverageInRegions.sh $tempFile /Users/acd13/Desktop/ATAC/Analysis/nucleosomeCalls/nucleosomeSignalInTFPeaks/${peakStage}TF100bpSummits/EEL3_H3DiffSignal
        rm ${tempFile}
	done
done
