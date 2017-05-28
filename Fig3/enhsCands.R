source('/Users/acd13/Desktop/ATAC/visualizing/gviz/gvizSupportFunctions_noLog.R')
setwd('/Users/acd13/Desktop/ATAC/visualizing/gviz/ourEnhCands')

nhr25Locus <- "X:13,007,000-13,016,000"
plotEE_v_L3_withDiffPeaks(nhr25Locus, comparisonDirection = "down", plotName="nhr25_EEL3Stages_25bpMeanEnrich_noLog")
plotL3_v_YA_withDiffPeaks(nhr25Locus, comparisonDirection = "up", plotName="nhr25_L3YAStages_25bpMeanEnrich_noLog")
plotAllStagesAndMetaPeaks(nhr25Locus, plotName="nhr25_allStages_25bpMeanEnrich_noLog")
plotAllStagesAndMetaPeaks(nhr25Locus, withChromStates=T, plotName="nhr25_allStages_25bpMeanEnrich_noLog_withChromHMM")

c54g6.3Locus <- "I:992295-1010294"
plotL3_v_YA_withDiffPeaks(c54g6.3Locus, comparisonDirection = "up", plotName="c54g6.3_L3YAStages_25bpMeanEnrich_noLog")
plotEE_v_L3_withDiffPeaks(c54g6.3Locus, comparisonDirection = "down", plotName="c54g6.3_EEL3Stages_25bpMeanEnrich_noLog")
plotAllStagesAndMetaPeaks(c54g6.3Locus, withChromStates=T, plotName="c54g6.3_allStages_25bpMeanEnrich_noLog_withChromHMM")
plotAllStagesAndMetaPeaks(c54g6.3Locus, plotName="c54g6.3_allStages_25bpMeanEnrich_noLog")

k12c11.3Locus <- "I:1,331,020-1,337,019"
plotL3_v_YA_withDiffPeaks(k12c11.3Locus, comparisonDirection = "down", plotName="k12c11.3_L3YAStages_25bpMeanEnrich_noLog")
plotEE_v_L3_withDiffPeaks(k12c11.3Locus, comparisonDirection = "up", plotName="k12c11.3_EEL3Stages_25bpMeanEnrich_noLog")
plotAllStagesAndMetaPeaks(k12c11.3Locus, withChromStates=T, plotName="k12c11.3_allStages_25bpMeanEnrich_noLog_withChromHMM")
plotAllStagesAndMetaPeaks(k12c11.3Locus, plotName="k12c11.3_allStages_25bpMeanEnrich_noLog")

tcc1 <- 'V:17,848,725-17,866,724'
#plotL3_v_YA_withDiffPeaks(tcc1, comparisonDirection = "up", plotName="tcc1_L3YAStages_25bpMeanEnrich_noLog", revStrand=TRUE)
plotEE_v_L3_withDiffPeaks(tcc1, comparisonDirection = "down", plotName="tcc1_EEL3Stages_25bpMeanEnrich_noLog", revStrand=TRUE)
plotAllStagesAndMetaPeaks(tcc1, withChromStates=T,plotName="tcc1_allStages_25bpMeanEnrich_noLog_withChromHMM", revStrand=TRUE)
plotAllStagesAndMetaPeaks(tcc1, plotName="tcc1_allStages_25bpMeanEnrich_noLog", revStrand=TRUE)

mlt8 <- "II:567,548-579,256"
plotEE_v_L3_withDiffPeaks(mlt8, comparisonDirection = "down", plotName="mlt8_EEL3Stages_25bpMeanEnrich_noLog", revStrand=TRUE)
plotL3_v_YA_withDiffPeaks(mlt8, comparisonDirection = "up", plotName="mlt8_L3YAStages_25bpMeanEnrich_noLog", revStrand=TRUE)
plotAllStagesAndMetaPeaks(mlt8, plotName="mlt8_allStages_25bpMeanEnrich_noLog", revStrand=TRUE)
plotAllStagesAndMetaPeaks(mlt8, withChromStates=T, plotName="mlt8_allStages_25bpMeanEnrich_noLog_withChromHMM", revStrand=TRUE)

swip10 <- "X:2,112,267-2,120,521"
plotEE_v_L3_withDiffPeaks(swip10, comparisonDirection = "down", plotName="swip10_EEL3Stages_25bpMeanEnrich_noLog", revStrand=TRUE)
plotL3_v_YA_withDiffPeaks(swip10, comparisonDirection = "up", plotName="swip10_L3YAStages_25bpMeanEnrich_noLog", revStrand=TRUE)
plotAllStagesAndMetaPeaks(swip10, plotName="swip10_allStages_25bpMeanEnrich_noLog", revStrand=TRUE)
plotAllStagesAndMetaPeaks(swip10, withChromStates=T, plotName="swip10_allStages_25bpMeanEnrich_noLog_withChromHMM", revStrand=TRUE)

gei13 <- "III:9,631,922-9,651,010"
plotEE_v_L3_withDiffPeaks(gei13, comparisonDirection = "down", plotName="gei13_EEL3Stages_25bpMeanEnrich_noLog", revStrand=TRUE)
plotL3_v_YA_withDiffPeaks(gei13, comparisonDirection = "up", plotName="gei13_L3YAStages_25bpMeanEnrich_noLog", revStrand=TRUE)
plotAllStagesAndMetaPeaks(gei13, plotName="gei13_allStages_25bpMeanEnrich_noLog", revStrand=TRUE)
plotAllStagesAndMetaPeaks(gei13, withChromStates=T, plotName="gei13_allStages_25bpMeanEnrich_noLog_withChromHMM", revStrand=TRUE)




