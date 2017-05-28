cd /Users/acd13/Desktop/ATAC/Analysis/differential/N2_dev/allButL1/differentialReports/combatCorrected/
# set up the files
tssFile="/Users/acd13/Desktop/geneDefinitions/refSeqGenes_TSSReplacedWithChenEtalTICs/tss.bed"
#parameters
# max dist is set below in a loop
# the FDR is hard coded below at 0.05
# there is no min FC

# I'll just loop through and do a number of things at once, then deal with them afterward
for MAX_DIST in 5000 10000 1e8 # max distance between tss and peak to be considered connected
do
    for stageComp in EE_vs_L3 EE_vs_YA L3_vs_YA
    do
        allPeaksDiffReport="allPeaks/${stageComp}_ComBatCorrected_EdgeR_q1.txt"
        # 1st 2 lines:
        #          Chr    Start      End       Conc     Conc_EE   Conc_L3          Fold      p-value          FDR          EEa          EEb          EEc          L3a          L3b        L3c
        # 546     chrI  2021187  2021551  5.7597847 -1.14106934  6.753736 -7.894805e+00 1.946486e-46 5.937562e-42 3.909441e-01 4.592282e-01 5.100978e-01     3.086133     4.669356  315.98591
        # I want Chr, Start, End, Conc_EE, Conc_L3, FDR
        # those are 2,3,4,6,7,10
        
        # step 1 is connecting every meta peak to a gene
        # luckily I've already done that for GO term enrichment, and all I have to do now is keep track of the peaks (previously I was only listing the genes)
        # this command comes from /Users/acd13/Desktop/ATAC/Analysis/differential/N2_dev/allButL1/differentialReports/getClosestRefChenTSS_narrowPeaks.sh
        # I've just tweaked it so that everything prints, instead of just the gene name
        
        connectFile="allPeaks/allPeaksConnectedToTssWithin${MAX_DIST}_${stageComp}Values_withStageConcentrations.bed"
        # this file has spaces in place of tabs, so get those to be single spaces only (currently any number), and then replace with tabs
        tail -n+2 ${allPeaksDiffReport} | tr -s ' ' | tr ' ' "\t" | cut -f 2,3,4,6,7,10 | closestBed -d -a - -b ${tssFile} | awk -v max="${MAX_DIST}" 'BEGIN{OFS="\t" } { if ( $NF <= max ) print $0 }' > ${connectFile}
        # an example line of the output
        # chrI	2021187	2021551	-1.14106934  6.753736	5.937562e-42	chrI	2026293	2026294	Y20F4.2	wormbase_tss	+	4743
        # which is peakChr, peakStart, PeakEnd, EE peak_concentration, L3_peak_concentration, FDR, TssChr, TssStart, TssEnd, GeneName, Class, Strand, distanceBTPeakAndTss
    done
done
