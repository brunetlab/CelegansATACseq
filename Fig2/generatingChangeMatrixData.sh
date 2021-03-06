cd /Users/acd13/Desktop/ATAC/Analysis/changingChromStates
#stage1="EE"
#stage2="L3"
stage1="L3"
stage2="YA"
outputDir="annotBedFile_${stage1}_vs_${stage2}"
mkdir -p ${outputDir}

# provide all peaks from diffBind
# split into 3 classes: up, down and neutral (should be able to get that from the peak labeling in predictAccessibility)
# for each class/file annotate with EE chromHMM state (just the big 7).
# then for each group of states , annotate with L3 chromHMM
# that should end up giving you 21 files with 7 lines (21 files=3 classes of peaks * 7 EE states, 7 lines = 7 L3 states)

# for ease I've just hardcoded the comparison file here because I think I'm only going to do this once or maybe 2x
peaks="/Users/acd13/Desktop/ATAC/Analysis/differential/N2_dev/allButL1/differentialReports/combatCorrected/allPeaks/${stage1}_vs_${stage2}_ComBatCorrected_EdgeR_q1.txt"
baseFilehandle="${outputDir}/${stage1}_v_${stage2}_fdr0.05"
# these files are a list of the chromHMM states (in bed format) of the main 7 states. One file path per line
firstChromStateFileList="/Users/acd13/Desktop/ATAC/Analysis/changingChromStates/${stage1}ChromHMMMergedStates.txt"
secondChromStateFileList="/Users/acd13/Desktop/ATAC/Analysis/changingChromStates/${stage2}ChromHMMMergedStates.txt"
# each of the above files is in the same order: Tx'd, TSS, act Enh, rep Enh, Rep'd, HetChrom, other
OVERLAP_PORTION=0.5

# the tail is to avoid the header
# I'll just write out three separate files, with 3 separate calls for ease
# first the no change peaks
# anything without a low enough FDR is considered no change
noChangePeaks="${baseFilehandle}_noChangePeaks.bed"
tail -n+2 ${peaks} | \
perl -lane 'if($F[9]>0.05){ print join("\t", ($F[1], $F[2], $F[3]));};' | \
sortBed -i - > ${noChangePeaks}

# now get down peaks
downPeaks="${baseFilehandle}_downPeaks.bed"
tail -n+2 ${peaks} | \
perl -lane 'if($F[9]<=0.05 && $F[7] < 0){ print join("\t", ($F[1], $F[2], $F[3]));};' | \
sortBed -i - > ${downPeaks}

# last the up peaks
upPeaks="${baseFilehandle}_upPeaks.bed"
tail -n+2 ${peaks} | \
perl -lane 'if($F[9]<=0.05 && $F[7] > 0){ print join("\t", ($F[1], $F[2], $F[3]));};' | \
sortBed -i - > ${upPeaks}

# now we'll do the same thing for each file, the comparison of them will come in an R script later
for t in noChangePeaks downPeaks upPeaks
do
    # it's stupid, but I'll recreate the names here, just so I can have that variable to change other names
    origFile="${baseFilehandle}_${t}.bed"
    # 
    while read firstLine ; do
		firstChromHMMSTATE=$(basename "${firstLine}" | sed 's/\.bed//g')
		firstStateFile="${baseFilehandle}_${t}_${firstChromHMMSTATE}.bed"
		coverageBed -a "${firstLine}" -b "${origFile}" | awk -v min="${OVERLAP_PORTION}" '{if ($NF > min) print $0}' | rev | cut -f 5- | rev > "${firstStateFile}"
		#
		while read secondLine ; do
			secondChromHMMSTATE=$(basename "${secondLine}" | sed 's/\.bed//g')
			finalFile="${baseFilehandle}_${t}_${firstChromHMMSTATE}_${secondChromHMMSTATE}.bed"
			coverageBed -a "${secondLine}" -b "${firstStateFile}" | awk -v min="${OVERLAP_PORTION}" '{if ($NF > min) print $0}' | rev | cut -f 5- | rev > "${finalFile}"
		done < "${secondChromStateFileList}"
	done < "${firstChromStateFileList}"    
done
