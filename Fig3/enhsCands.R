source('/Users/acd13/Desktop/ATAC/visualizing/gviz/gvizSupportFunctions_noLog.R')
setwd('/Users/acd13/Desktop/ATAC/visualizing/gviz/ourEnhCands')

nhr25Locus <- "X:13,007,000-13,016,000"
plotAllStagesAndMetaPeaks(nhr25Locus, plotName="nhr25_allStages_25bpMeanEnrich_noLog")

c54g6.3Locus <- "I:992295-1010294"
plotAllStagesAndMetaPeaks(c54g6.3Locus, plotName="c54g6.3_allStages_25bpMeanEnrich_noLog")

mlt8 <- "II:567,548-579,256"
plotAllStagesAndMetaPeaks(mlt8, plotName="mlt8_allStages_25bpMeanEnrich_noLog", revStrand=TRUE)

swip10 <- "X:2,112,267-2,120,521"
plotAllStagesAndMetaPeaks(swip10, plotName="swip10_allStages_25bpMeanEnrich_noLog", revStrand=TRUE)

gei13 <- "III:9,631,922-9,651,010"
plotAllStagesAndMetaPeaks(gei13, plotName="gei13_allStages_25bpMeanEnrich_noLog", revStrand=TRUE)
