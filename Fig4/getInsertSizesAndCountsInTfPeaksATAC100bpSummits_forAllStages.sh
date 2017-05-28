peakStage=$1
controlAtacInsertSiteFile="/Users/acd13/Desktop/ATAC/insertSites/N2dev/N2_dev_gDNA.raw.srt.filt.nodup.srt.shifted.insertSites.bed"
cd "/Users/acd13/Desktop/ATAC/Analysis/tfs/qualityTfs/insertSizeStatsAtPeaks/ATACSignalSummitsInTFs"
mkdir -p ${peakStage}
cd ${peakStage}

perPeakDirSize="rawResults/perTfPeak/insertSizes"
perPeakDirCounts="rawResults/perTfPeak/counts"

for dir in rawResults ${perPeakDirAllCounts} ${perPeakDirSize}
do
	mkdir -p ${dir}
done

for stage in EE L3 YA
do
	mkdir -p "${perPeakDirSize}/${stage}" 
	mkdir -p "${perPeakDirCounts}/${stage}" 
	# we'll also need the total number of reads for the stage specific and control counts, to make enrichment calculations
	# but we already did that, so we'll take advantage of this loop to  link to those
	ln -s /Users/acd13/Desktop/ATAC/Analysis/tfs/qualityTfs/insertSizeStatsAtPeaks/${peakStage}/rawResults/${stage}AtacInsertSiteTotal.txt ./rawResults/${stage}AtacInsertSiteTotal.txt
done
# and for the gDNA
ln -s /Users/acd13/Desktop/ATAC/Analysis/tfs/qualityTfs/insertSizeStatsAtPeaks/${peakStage}/rawResults/gDNAControlAtacInsertSiteTotal.txt ./rawResults/gDNAControlAtacInsertSiteTotal.txt

mkdir -p "${perPeakDirCounts}/gDNAControl" 

function getATACData {
    tfPeakFile=${1}
    perPeakDirSize=${2}
    perPeakDirCounts=${3}
    controlAtacInsertSiteFile=${4}
    prefix=$(basename ${tfPeakFile} | perl -lane '@nameParts = split/\.narrowPeak\.gz/, $F[0]; print $nameParts[0];')
    echo "working on ${prefix}"
    tempFile="${prefix}.tmp"
    # b/c these are narrowPeaks, they have the summits as the last column
    sortBed -i ${tfPeakFile} | slopBed -b 45 -i - -g $CE10_GSIZE > ${tempFile} # this results in 100bp windows b/c the summits are 10bp windows
    
     for stage in EE L3 YA
     do
        stageSpecificAtacInsertSiteFile="/Users/acd13/Desktop/ATAC/insertSites/N2dev/${stage}_allReps.adjusted.insertSites.bed.gz"
        intersectBed -wb -a ${stageSpecificAtacInsertSiteFile} -b ${tempFile} | \
        perl -lane '$val=abs($F[4]); print( join("\t", ($F[6], $F[7], $F[8], $val)));' | \
        sortBed -i - | bedtools groupby -i - -g 1,2,3 -c 4 -o median | \
        perl -lane 'print(join("_", ($F[0], $F[1], $F[2]))."\t".$F[3]);' > "${perPeakDirSize}/${stage}/${prefix}.txt"
        
        coverageBed -a ${stageSpecificAtacInsertSiteFile} -b ${tempFile} | sortBed -i - | \
        perl -lane 'print(join("_", ($F[0], $F[1], $F[2]))."\t".$F[3]);' > "${perPeakDirCounts}/${stage}/${prefix}.narrowPeak"
    done
    coverageBed -a ${controlAtacInsertSiteFile} -b ${tempFile} | sortBed -i - | \
    perl -lane 'print(join("_", ($F[0], $F[1], $F[2]))."\t".$F[3]);' > "${perPeakDirCounts}/gDNAControl/${prefix}.narrowPeak"
    
    rm ${tempFile}
}
export -f getATACData

parallel 'getATACData {}' ::: /Users/acd13/Desktop/ATAC/Analysis/tfs/qualityTfs/peaksToUse/${peakStage}/*.narrowPeak.gz_10bp${peakStage}ATACSummit.bed.gz ::: ${perPeakDirSize} ::: ${perPeakDirCounts} ::: ${controlAtacInsertSiteFile}

Rscript /Users/acd13/Desktop/ATAC/Analysis/tfs/plottingAtacInsertSizesAndEnrichmentInTfPeakATAC100bpSummits_withChanges.R ${peakStage}
