ln -s /Users/acd13/Desktop/ATAC/macs2Peaks/usingInserts_75bpShift/N2_dev/withAllReps/allRepsAllButL1MetaPeaks_smartMerged300bpSummitDist.bed ./consensus_peaks.bed

# The next 2 lines were repeated for EE and L3 as well
bash /Users/acd13/Desktop/ATAC/Analysis/Enrichments/myChromHMMPeakAnnotater.sh /Users/acd13/Desktop/ATAC/macs2Peaks/usingInserts_75bpShift/N2_dev/withAllReps/allRepsAllButL1MetaPeaks_smartMerged300bpSummitDist.bed YA
for f in *.bed; do name=$(echo $f | cut -f 4- -d'_'| cut -f 1 -d'.'); while read line; do echo -e "${line}\t${name}"; done < ${f} ; done | sortBed -i - > YAChromHMMAnnotated_consensusPeaks.bed

bash /Users/acd13/Dropbox/Scripts/MostUsed/geneAnnotaters/myPeakAnnotaterRefseqWithChenTSSNoMirs.sh consensus_peaks.bed
for f in *.bed; do name=$(echo $f | cut -f3- -d'_'| rev | cut -f2- -d'.'|rev); while read line; do echo -e "${line}\t${name}"; done < ${f} ; done | sortBed -i - > gene_location_annotated_consensusPeaks.bed

# same was done for L3 v YA
cat /Users/acd13/Desktop/ATAC/Analysis/differential/N2_dev/allButL1/differentialReports/combatCorrected/allPeaks/EE_vs_L3_ComBatCorrected_EdgeR_q1.tsv | cut -f 2-4,8,10 | tail -n+2 | intersectBed -wa -wb -a - -b consensus_peaks.bed | sortBed -i - > differential_annots/EE_vs_L3_ComBatCorrected_EdgeR_FC_FDR.bed
# I added the intersection b/c DiffBind merged overlapping peaks

closestBed -D ref -a consensus_peaks.bed -b ~/Desktop/geneDefinitions/refSeqGenes_TSSReplacedWithChenEtalTICs/tss.bed | cut -f 1-4,8,11| sort -u | sortBed -i - > nearestTss_consensusPeaks.bed

python combiner.py nearestTss_consensusPeaks.bed EEChromHMMStateAnnotatedPeaks/EEChromHMMAnnotated_consensusPeaks.bed L3ChromHMMStateAnnotatedPeaks/L3ChromHMMAnnotated_consensusPeaks.bed YAChromHMMStateAnnotatedPeaks/YAChromHMMAnnotated_consensusPeaks.bed annotated_peaks/gene_location_annotated_consensusPeaks.bed differential_annots/EE_vs_L3_ComBatCorrected_EdgeR_FC_FDR.bed differential_annots/L3_vs_YA_ComBatCorrected_EdgeR_FC_FDR.bed
head -1 fullyAnnotated_consensusAtacPeaks.bed > temp
tail -n+2 fullyAnnotated_consensusAtacPeaks.bed | sortBed -i - >> temp
mv temp fullyAnnotated_consensusAtacPeaks.bed
