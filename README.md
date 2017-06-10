# CelegansATACseq
The code for Daugherty, et al 2017 - Chromatin accessibility dynamics reveal novel functional enhancers in C. elegans


Figure|Panel|Comments|Program name
---|---|---|---
1|A/B|General pipeline|ATACPipeline_combined.sh
|||pyadapter_trim.py
|||random_split_fastq.pl
---|---|---|---
	|	|get consensus peaks|reRunningWithoutL1s.sh
---|---|---|---
	|	|	|splitMACS2SubPeaks.pl
---|---|---|---
	|	|	|smart_merge.py
---|---|---|---
	|C|Differential, as well as 1B|30JanDiffBind_allButL1_metaPeaks_withGDNA.R
---|---|---|---
	|D & E|Gene plots|examples.R
---|---|---|---
	|	|support for the gene plots|gvizSupportFunctions_noLog.R
---|---|---|---
	|F & G|How genes were connected to peaks; genes were then sorted using command line sort, and gene lists copied over to Gorilla|connectingPeaksAndTss_withConcentrations.sh
---|---|---|---
	|	|Many of these combinations were ad hoc, so this is harder to understand|selectTerms_ggPlot2_EEvL3_allConnectedATACPeakTotalChange_GOProcess_splitup_selectGroupedTerms.R
---|---|---|---
|||
---|---|---|---			
S1|A|Plots from Picard in the General pipeline|From above (1 A/B)
---|---|---|---
	|B||jaccardOfRegions.R
---|---|---|---
	|C/D|Volcano plots|plotDifferences.24Mar2015_volcano.R
---|---|---|---
	|E&F|Many of these combinations were ad hoc, so this is harder to understand|selectTerms_ggPlot2_L3vYA_allConnectedATACPeakTotalChange_GOProcess_splitup_selectGroupedTerms
---|---|---|---
	|	|	|
---|---|---|---
S2|A/B|description and how each of the remaining scripts was called; mini pipeline|whatWasDone
---|---|---|---
	|	|called in the above|trim_galore_SE.sh
---|---|---|---
	|	|called in the above|mapqFilterAndStrandSpecificSplit.sh
---|---|---|---
	|	|called in the above|GROseq_alignment_bowtie2.sh
---|---|---|---
	|	|called in the above|getConsensusPeaks_2BioReps_version4.sh
---|---|---|---
	|	|called in the above|deDup.sh
---|---|---|---
	|	|Plotting of the ATAC and GRO-seq data|plottingPromAccessibilityVsGroSignalRefChen.R
---|---|---|---
	|	|	|
---|---|---|---
	|C/D|All the scripts for RNA-seq analysis|runningUcscAlignments.sh
---|---|---|---
	|	|All the scripts for RNA-seq analysis|topHat2_RNAseq_alignmentForModEncodeData_36bp_PE.sh
---|---|---|---
	|	|All the scripts for RNA-seq analysis|topHat2_RNAseq_alignmentForModEncodeData_76bpReads_SE.sh
---|---|---|---
	|	|All the scripts for RNA-seq analysis|quantifyTophat2ResultsWithHtSeq_RefSeq.sh
---|---|---|---
	|	|Plotting of the ATAC and RNA-seq data|See plottingPromAccessibilityVsGroSignalRefChen.R
---|---|---|---
	|	|	|
---|---|---|---
2|A|Step by step of what was run|ExactlyWhatIdid.txt
---|---|---|---
	|	|all of the various scripts|chromHMM_scripts
---|---|---|---
	|	|Generating the null distribution|getEnrichments_justBash.sh
---|---|---|---
	|	|Measuring intersection with the null|createShuffledBedDirectory_mappableSubtractBlacklistAndGDNAPeaks_parallel.sh
---|---|---|---
	|	|actually generating the bar plots|plottingAllTogetherInTheirStates.R
---|---|---|---
	|B|Standard NGS plot code|Example command: ngs.plot.r -G ce10 -R bed -C ConfigVsInput.txt -O L3_HM_SignalVsInput_atL3ATACPeaks_localScale -P 0 -IN 1 -GO max -N 0.33 -VLN 0 -CO blue:red
---|---|---|---
	|	|How the Chen TSSs were used to replace the canonical TSS|replaceFeaturesWherePossible_byName.py
---|---|---|---
	|C/D|Getting the data|generatingChangeMatrixData.sh
---|---|---|---
	|	|plotting|chreatingChromHMMStateChangeHeatMaps.R
---|---|---|---
	|E|Example plot|examples.R
---|---|---|---
	|	|	|
---|---|---|---
S3|A-C|See ChromHMM code above|N/A
---|---|---|---
	|	|	|plottingStackedBarPlots.R
---|---|---|---
	|D/E|Getting the data|splitAndCombineHiHMMStates.sh
---|---|---|---
	|	|	|compare_hiHMMToMyChromHMM.sh
---|---|---|---
	|	|	|runningCompares.sh
---|---|---|---
	|	|plotting|jaccardOfChromatinStatePreds.R
---|---|---|---
	|F|Standard NGS plot code|See above
---|---|---|---
	|G/H|See 2C/D|N/A
---|---|---|---
	|	|	|
---|---|---|---
S4|A|See 2A
---|---|---|---
	|B|Generating Null distribution|generatingAllDistalNegsOnScg3_withoutProms.sh
---|---|---|---
	|	|	|plottingConservationScoresWithNulls.R
---|---|---|---
	|C|	|plottingConservationScoresVsHistoneMods_noNulls_medians.R
---|---|---|---
	|D|See examples.R|for general approach	
---|---|---|---
	|	|	|
---|---|---|---
3/S5|all plots|	|enhsCands.R
---|---|---|---
	|	|	|
---|---|---|---
4|A|finding known motifs|findMotifsGenome_sizeGiven_narrowPeak_noL1sMetaPeakBackground_novelMotifsLO9.sh
---|---|---|---
	|B|Labeling peaks|makeEEL3Labels.sh
---|---|---|---
	|	|Getting the counts of each motif in the peaks|getMotifCounts_EEvL3DeNovoToo.sh
---|---|---|---
	|	|Final prep of the input data|condenseCountsPerPeaks_withPreDoneNames.py
---|---|---|---
	|	|building then predicting from the model|predictingAccessibilityChangesWithHughesAndDeNovoMotifs_gbmMetricBalAcc_LO9Peaks_loopingOverInteractionDepth_noTSSDist.R
---|---|---|---
	|C|This isn't the exact code used, but this is a functio used in all predictors, so it's identical to what was then modified to make prettier for plotting|plotting_rel_importance.R
---|---|---|---
	|D|Getting the data|getInsertSizesAndCountsInTfPeaksATAC100bpSummits_forAllStages.sh
---|---|---|---
	|	|Processing the data|processingAtacInsertSizesAndEnrichmentInTfPeak100bpSummits.R
---|---|---|---
	|	|supporting the above|ProcessingAtacInsertSizesAndEnrichmentInTfPeaks_supportFunctions.r
---|---|---|---
	|	|Plotting|plottingATACSignalAsFunctionOfATACFragmentSize.R
---|---|---|---
	|E|just the plotting|plotting_eor1_insertSize_hist.R
---|---|---|---
	|F|wrapper script for DANPOS|runningDanpos.sh
---|---|---|---
	|	|called above|wigToBedGraph.py
---|---|---|---
	|	|Calculate the differences|nucDiffSignalInTF100bpSummits.sh
---|---|---|---
	|	|called above|getH3NucleosomeEEL3DiffCoverageInRegions.sh
---|---|---|---
	|	|plotting|plottingInferredNucEEL3DiffSignalInTF100bpSummits.R
---|---|---|---
	|	|called above|plottingAtacInsertSizesAndEnrichmentInTfPeaks_supportFunctions.R
---|---|---|---
	|	|	|
---|---|---|---
S6|A|See 4A for example	|
---|---|---|---
	|B|See 4A for example	|
---|---|---|---
	|C|de  novo motifs were called using standard homer settings. See methods for details|
---|---|---|---
	|D|plotting accuracy of model|plottingBalancedAccuracy.R
---|---|---|---
	|	|	|
---|---|---|---
S7|A|See S6C for example|	
---|---|---|---
	|B|See 2A for example of approach|
---|---|---|---
	|C/D|Plotting boxplots for each factor|plottingATACSignalAndFragmentSizeInTFChIPPeaks_v2.R
---|---|---|---
	|E|See 4E	|
---|---|---|---
	|F|Generating the data|bedtoolsFisherTestOneVsManyWrapper.sh
---|---|---|---
	|	|Plotting|plot_fishers_test_ORs_L3inL3.R
---|---|---|---
	|	|	|
---|---|---|---
Sup Table 3|	|Mostly copy and paste of how I pieced together the annotation files|annotating_consensus_atac_peaks.sh
---|---|---|---
	|	|used in the above to rename some things|combiner.py
---|---|---|---
