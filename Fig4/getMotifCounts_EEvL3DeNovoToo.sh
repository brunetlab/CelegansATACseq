# this is similar to the original getMotifCounts, but I've added de novo motifs discovered by homer in EE vs L3 differential atac peaks, and mapped by Homer


# this takes mapped PWMs from Homer, in bed format, (e.g. output of scanMotifGenomeWide.pl) and a set of regions
# The bed format from Homer has the motif name in column 4, and the logOdds in column 5, both 1 based
# This pulls out the motifs instances that entirely intersect the peaks, and summarizes the number of instances of every motif found in each peak
# The for loop does this at a range of log odds to allow you to play around with different stringencies of motif mapping

loMin=9
denovoUp='/Users/acd13/Desktop/ATAC/Analysis/Motifs/ATAC_differential/N2_dev/combatCorrectedFdr0.05/peaks/motifs/metaPeakBackground/EE_vs_L3_ComBatCorrected_EdgeR_q0.05_up/homerResults/allHomerMotifs.mappedToCe10.bed.gz'
denovoDown='/Users/acd13/Desktop/ATAC/Analysis/Motifs/ATAC_differential/N2_dev/combatCorrectedFdr0.05/peaks/motifs/metaPeakBackground/EE_vs_L3_ComBatCorrected_EdgeR_q0.05_down/homerResults/allHomerMotifs.mappedToCe10.bed.gz'
mappedHughesMotifs="hughesCePWMSForHomerLO5.motif_genomewide_ce10.srt.bed.gz"
mappedMotifs="hughesCePWMSForHomerLO5_plusUpAndDownHomerResultsFrom_EE_v_L3_differentialATACPeaks_genomewide_ce10.srt.bed.gz"
peaks="/Users/acd13/Desktop/ATAC/Analysis/differential/N2_dev/allButL1/differentialReports/combatCorrected/allPeaks/allPeaksConsidered.bed"
outputFileHandle="motifCounts/hughesMotifAndDeNovoCountsInMetaAtacPeaks"

cd /Users/acd13/Desktop/ATAC/Analysis/predictingAccessChangeWithMotifs

# now add in de novo motifs and sort everything
#zcat "${mappedHughesMotifs}" "${denovoUp}" "${denovoDown}" | sort -k 1,1 -k2,2n | gzip -c > "${mappedMotifs}"
# probably from homer mapping, I ran into one or 2 lines that didn't have intelligble bed intervals, so I'll just filter those
zcat "${mappedMotifs}" | perl -lane 'print unless ($F[2] - $F[1] < 2);' | gzip -c > temp
mv temp "${mappedMotifs}"

# I was having issues with a single motif being listed 2x, so I added the sort. Effectively what I do is sort by peak position, and then by motif name, to ensure that all instances of the motif within a peak are next to each other
intersectBed -f 1 -wb -a ${mappedMotifs} -b $peaks | awk -v thresh="${loMin}" '{if ($5>thresh){print}}' | sort -k 7,7 -k8,8n -k4,4 | groupBy -i - -g 7,8,9,4 -c 4 -o count > "${outputFileHandle}_LO${loMin}.bed"
python condenseCountsPerPeaks.py -f "${outputFileHandle}_LO${loMin}.bed" > "${outputFileHandle}_LO${loMin}.motifOccurenceList.txt"
