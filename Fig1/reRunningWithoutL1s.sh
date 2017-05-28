# Call raw peaks from the merged reads
cd /Users/acd13/Desktop/ATAC/insertSites/N2dev/shiftedForMacs2/chromEndsRemoved
macs2 callpeak -t /Users/acd13/Desktop/ATAC/insertSites/N2dev/shiftedForMacs2/chromEndsRemoved/allRepsAllButL1.adjusted.insertSites.75bpShiftedForMACS.bed -f BED -n /Users/acd13/Desktop/ATAC/macs2Peaks/usingInserts_75bpShift/N2_dev/withAllReps/allRepsAllButL1.MACS2.1_q5e-2 -g 9e7 -p 5e-2 --nomodel --extsize 75 -B --keep-dup all --call-summits
# move things around
cd /Users/acd13/Desktop/ATAC/macs2Peaks/usingInserts_75bpShift/N2_dev/withAllReps
mkdir allRepsAllButL1
for file in  allRepsAllButL1.MACS2.1_q5e-2* ; do mv $file allRepsAllButL1/; done
cd allRepsAllButL1
# Split sub-peaks, remove blacklists, and gDNA intersecting peaks
perl /Users/acd13/Dropbox/Scripts/ATAC_scripts/splitMACS2SubPeaks.pl allRepsAllButL1.MACS2.1_q5e-2_peaks.narrowPeak > allRepsAllButL1.MACS2.1_q5e-2_peaks.splitSubsPeaks.narrowPeak  
intersectBed -v -a allRepsAllButL1.MACS2.1_q5e-2_peaks.splitSubsPeaks.narrowPeak -b /Volumes/extra/Genomic_files/worm/ce10/ce10-blacklist.bed > allRepsAllButL1.MACS2.1_q5e-2_peaks.splitSubsPeaks.blacklistCleared.narrowPeak
intersectBed -u -a allRepsAllButL1.MACS2.1_q5e-2_peaks.splitSubsPeaks.narrowPeak -b /Volumes/extra/Genomic_files/worm/ce10/ce10-blacklist.bed > allRepsAllButL1.MACS2.1_q5e-2_peaks.splitSubsPeaks.intersectingBlacklist.narrowPeak
intersectBed -f 0.2 -v -a allRepsAllButL1.MACS2.1_q5e-2_peaks.splitSubsPeaks.blacklistCleared.narrowPeak -b /Users/acd13/Desktop/ATAC/macs2Peaks/usingInserts_75bpShift/N2_dev/withAllReps/gDNAMasked_MACS5e-2_0.5Overlap/N2_gDNA.filt.nodup.adjusted.macs2.1q5e-2_peaks.narrowPeak > allRepsAllButL1.MACS2.1_q5e-2_peaks.splitSubsPeaks.blacklistCleared.0.2OverlapGDNAMasked.narrowPeak

# All of the above was from the initial attempt that didn't use the psuedo reps, 
# That was combine with the command below to get the actual consensus peaks
bash /Users/acd13/Desktop/ATAC/macs2Peaks/usingInserts_75bpShift/N2_dev/withAllReps/callMetaUnion_L1s.sh
# The final file I cared about from that is:
# allRepsAllButL1.MACS2.1_q5e-2_splitSubsPeaks_blacklistCleared_consensusPeaks_gDNA0.2PeaksRemoved.narrowPeak
# get the summits of these remaining peaks for use with the smartMerge from Daniel Kim
intersectBed -f 1 -u -a allRepsAllButL1.MACS2.1_q5e-2_summits.bed -b allRepsAllButL1.MACS2.1_q5e-2_splitSubsPeaks_blacklistCleared_consensusPeaks_gDNA0.2PeaksRemoved.narrowPeak > allRepsAllButL1.MACS2.1_q5e-2_peaks_splitSubsPeaks_blacklistCleared_consensusPeaks_gDNA0.2PeaksRemoved.summits.bed
#
# create summary file for smart merge
echo -e "allRepsAllButL1.MACS2.1_q5e-2_splitSubsPeaks_blacklistCleared_consensusPeaks_gDNA0.2PeaksRemoved.narrowPeak\t/Users/acd13/Desktop/ATAC/macs2Peaks/usingInserts_75bpShift/N2_dev/withAllReps/allRepsAllButL1/allRepsAllButL1.MACS2.1_q5e-2_peaks_splitSubsPeaks_blacklistCleared_consensusPeaks_gDNA0.2PeaksRemoved.summits.bed" > forSmartMerge_withoutL1s.tsv
echo -e "/Users/acd13/Desktop/ATAC/macs2Peaks/usingInserts_75bpShift/N2_dev/withAllReps/gDNAMasked_MACS5e-2_0.5Overlap/EE_75bpMACS2Shift_mergedReps.MACS2.1_q5e-2_splitSubsPeaks_consensusPeaks_gDNA0.2PeaksRemoved.narrowPeak\t/Users/acd13/Desktop/ATAC/macs2Peaks/usingInserts_75bpShift/N2_dev/withAllReps/gDNAMasked_MACS5e-2_0.5Overlap/summits/EE_75bpMACS2Shift_mergedReps.MACS2.1_q5e-2_splitSubsPeaks_consensusPeaks_gDNA0.2PeaksRemoved.summits.bed" >>forSmartMerge_withoutL1s.tsv
echo -e "/Users/acd13/Desktop/ATAC/macs2Peaks/usingInserts_75bpShift/N2_dev/withAllReps/gDNAMasked_MACS5e-2_0.5Overlap/L3_75bpMACS2Shift_mergedReps.MACS2.1_q5e-2_splitSubsPeaks_consensusPeaks_gDNA0.2PeaksRemoved.narrowPeak	/Users/acd13/Desktop/ATAC/macs2Peaks/usingInserts_75bpShift/N2_dev/withAllReps/gDNAMasked_MACS5e-2_0.5Overlap/summits/L3_75bpMACS2Shift_mergedReps.MACS2.1_q5e-2_splitSubsPeaks_consensusPeaks_gDNA0.2PeaksRemoved.summits.bed" >>forSmartMerge_withoutL1s.tsv
echo -e "/Users/acd13/Desktop/ATAC/macs2Peaks/usingInserts_75bpShift/N2_dev/withAllReps/gDNAMasked_MACS5e-2_0.5Overlap/YA_75bpMACS2Shift_mergedReps.MACS2.1_q5e-2_splitSubsPeaks_consensusPeaks_gDNA0.2PeaksRemoved.narrowPeak	/Users/acd13/Desktop/ATAC/macs2Peaks/usingInserts_75bpShift/N2_dev/withAllReps/gDNAMasked_MACS5e-2_0.5Overlap/summits/YA_75bpMACS2Shift_mergedReps.MACS2.1_q5e-2_splitSubsPeaks_consensusPeaks_gDNA0.2PeaksRemoved.summits.bed" >>forSmartMerge_withoutL1s.tsv
# And finally run the smart merge
python /Users/acd13/Dropbox/Scripts/ATAC_scripts/smart_merge.py forSmartMerge_withoutL1s.tsv 150 /Users/acd13/Desktop/ATAC/macs2Peaks/usingInserts_75bpShift/N2_dev/withAllReps/allRepsAllButL1MetaPeaks_smartMerged150bpSummitDist.bed
mv forSmartMerge_withoutL1s.tsv /Users/acd13/Desktop/ATAC/macs2Peaks/usingInserts_75bpShift/N2_dev/withAllReps/
