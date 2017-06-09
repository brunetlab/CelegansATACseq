source('/Users/acd13/Desktop/ATAC/visualizing/gviz/gvizSupportFunctions_noLog.R')
setwd('/Users/acd13/Desktop/ATAC/visualizing/gviz/differentialExamples/')


cav1 <- 'IV:9,770,419-9,774,512'
plotAllStagesAndMetaPeaks(cav1, 10, plotName="cav1_allStages_10bpMeanEnrich_noLog")

daf12 <- 'X:10,629,613-10,667,976'
plotAllStagesAndMetaPeaks(daf12, 25, plotName="daf12_allStages_25bpMeanEnrich_noLog")


kin1 <- 'I:14,916,548-14,950,960'
plotEE_v_L3_withDiffPeaks(kin1, 25, 'down', withChromStates=TRUE, plotName="kin1_EEvL3_25bpMeanEnrich_withChromHMM_noLog")
