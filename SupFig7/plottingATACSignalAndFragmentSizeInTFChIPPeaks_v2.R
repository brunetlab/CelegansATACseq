stage <- 'L3'
stageInd <- 2
color <- 'goldenrod2'
wd <- paste0('/Users/acd13/Desktop/ATAC/Analysis/tfs/qualityTfs/insertSizeStatsAtPeaks/TFPeakSummits/',stage)

plottingTFStats <- function(toPlot, plotName, stage, color){
  pdf(paste0(stage, '_ATACSeq', plotName, 'In', stage, 'TFPeaks.pdf'), height = 5, width = 10)
    ## Increase bottom margin to make room for rotated labels
    par(mar = c(7, 4, 4, 2) + 0.1)
    ## Create plot with no x axis and no x axis label
    boxplot(toPlot
            , col=color
            , notch=T
            , outline=F
            , ylab=paste0('ATAC-seq ',  plotName, ' in TF ChIP-seq peaks')
            , xlab=""
            , xaxt = "n"
            , las = 1
    )
    
    ## Set up x axis with tick marks alone
    axis(1, at = 1:length(toPlot), labels = FALSE)
    ## Plot x axis labels at default tick marks
    text(1:length(toPlot)
         , par("usr")[3] - 0.5
         , srt = 45
         , adj = 1
         , labels = names(toPlot)
         , xpd = TRUE
    )
    ## Plot x axis label at line 6 (of 7)
    mtext(1, text = paste0(stage, ' TF ChIP-seq peaks'), line = 6)
    # When plotting the x axis labels, we use srt = 45 for text rotation angle, 
    # adj = 1 to place the right end of text at the tick marks, 
    # and xpd = TRUE to allow for text outside the plot region. 
    # You can adjust the value of the 0.25 offset as required to move the 
    # axis labels up or down relative to the x axis
  dev.off()  
}



setwd(wd)
size <- readRDS('perTFPeakAllStageATACInsertSize.list.rds')
signal <- readRDS('perTFPeakAllStageATACEnrichment.list.rds')

stageSpecificSize <- size[[stageInd]]
### 19 Feb 2017 I'm removing the DCC EOR-1 b/c it is not a trusted data source
stageSpecificSize[grep('DCC3160', names(stageSpecificSize))] <- NULL
names(stageSpecificSize) <- unlist(lapply(strsplit(names(stageSpecificSize), "_"), `[[`, 2))
sizeOrderedPlotList <- stageSpecificSize[rev(order(
                                                unlist(lapply(stageSpecificSize, function(x){median(unlist(x))})),
                                                unlist(lapply(stageSpecificSize, function(x){quantile(unlist(x), 0.75)}))
                                                ))]
stageSpecificSignal <- signal[[stageInd]]
stageSpecificSignal[grep('DCC3160', names(stageSpecificSignal))] <- NULL
names(stageSpecificSignal) <- unlist(lapply(strsplit(names(stageSpecificSignal), "_"), `[[`, 2))
signalOrderedPlotList <- stageSpecificSignal[order(
  unlist(lapply(stageSpecificSignal, function(x){median(unlist(x))})),
  unlist(lapply(stageSpecificSignal, function(x){quantile(unlist(x), 0.75)}))
  ,decreasing = T
  )
  ]


plottingTFStats(sizeOrderedPlotList, "revOrderedSize", stage, color)
plottingTFStats(signalOrderedPlotList, "orderedSignal", stage, color)


