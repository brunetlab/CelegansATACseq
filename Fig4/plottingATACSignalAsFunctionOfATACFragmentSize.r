stage <- 'L3'
stageInd <- 2
library('scales')
setwd(paste0('/Users/acd13/Desktop/ATAC/Analysis/tfs/qualityTfs/insertSizeStatsAtPeaks/TFPeakSummits/',stage))
size <- readRDS('perTFPeakAllStageATACInsertSize.list.rds')
signal <- readRDS('perTFPeakAllStageATACEnrichment.list.rds')

stageSpecificSsize <- size[[stageInd]]
stageSpecificSignal <- signal[[stageInd]]

medianSizes <- unlist(lapply(stageSpecificSsize, median))
medianSignal <- unlist(lapply(stageSpecificSignal, median))

if(all(names(medianSizes) != names(medianSignal))){
  warning("Size and signal aren't in the same order")
}

eor1Indicies <- grep('EOR-1', names(medianSizes))

pdf(paste0('ATACSeqSignalAsFunctionOfATACFragmentSizeIn', stage, 'TFPeaks.pdf'), height = 5, width = 5)
  plot(medianSizes, medianSignal
     , pch=16
     , col=alpha('goldenrod2', 0.7)
     , cex=1.5
     , ylim=c(0,max(medianSignal))
     , xlim=c(100,max(medianSizes)*1.1)
     , xlab='Median ATAC-seq fragment size in TF ChIP-seq peaks'
     , ylab='Median ATAC-seq signal in TF ChIP-seq peaks')
points(medianSizes[eor1Indicies], medianSignal[eor1Indicies], pch=1, col = 'black', cex=1.5)
legend('topright', legend = c("EOR-1", "Other L3 factors"), col = c('black', 'grey'), pch=c(1,1,16))
dev.off()  




