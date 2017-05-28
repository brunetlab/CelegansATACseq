tail -n+2 /Users/acd13/Desktop/ATAC/Analysis/differential/N2_dev/allButL1/differentialReports/combatCorrected/allPeaks/EE_vs_L3_ComBatCorrected_EdgeR_q1.txt | \
perl -lane 'if($F[9]>0.05){$label=0}elsif($F[7]<0){$label=-1}else{$label=1}; print join("\t", ($F[1], $F[2], $F[3], join("_", ($F[1], $F[2], $F[3])), $label));' | \
sortBed -i - | cut -f4- > metaAtacPeaks_EE_L3_fdr0.05Labels.txt

