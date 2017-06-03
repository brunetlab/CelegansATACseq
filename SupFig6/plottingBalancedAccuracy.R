# read in data previously generated
setwd('/Users/acd13/Desktop/ATAC/Analysis/predictingAccessChangeWithMotifs/usingGBM/predicting_metaAtacPeaks_EE_L3_fdr0.05Labels_withLO9HughesAndDeNovoMotifs_peaks_noTssDist')

trainingStats <- 'predicting_metaAtacPeaks_EE_L3_fdr0.05Labels_withLO9HughesAndDeNovoMotifs_peaks_interactionDepth_8_metricBalAcc_trainSetConfusionMatrixAndStats.txt'
testingStats <- 'predicting_metaAtacPeaks_EE_L3_fdr0.05Labels_withLO9HughesAndDeNovoMotifs_peaks_interactionDepth_8_metricBalAcc_testSetConfusionMatrixAndStats.txt'

trainingData <-  read.table(trainingStats, header=F, skip=10, row.names=2, stringsAsFactors=F)[,-1]
testingData <-  read.table(testingStats, header=F, skip=10, row.names=2, stringsAsFactors=F)[,-1]

# unfortunately I just have to add the column names, but it was copied from here
# http://www.inside-r.org/node/86995 the help page for the function used to generate this data
statsColumnNames <- c("sensitivity", "specificity", "positive predictive value", "negative predictive value"
                      , "prevalence", "dection rate", "detection prevalence", "balanced accuracy")

colnames(trainingData) <- colnames(testingData) <- statsColumnNames

pdf('testingSet_predictionStats.pdf', width=4, height=5)
  barplot(as.matrix(t(testingData[c(3,2,1),"balanced accuracy"]))
                  ,beside=T
                  , ylab='Balanced accuracy'
                  , ylim=c(0,1)
                  , xlab='ATAC peak dynamics'
                  ,col = c('darkorchid4','grey', 'goldenrod2')
                  ,names.arg = c("Embryo", "No change", "Larval")
                    )
abline(h=0.5, col='red')

barplot(as.matrix(t(testingData[c(3,2,1),c("sensitivity", "specificity")]))
                ,beside=T
                , ylab='Balanced accuracy'
                , ylim=c(0,1)
                , xlab='ATAC peak dynamics'
                ,col = c('darkorange1', 'cyan4')
                ,names.arg = c("Embryo", "No change", "Larval")
)
legend('top', legend=c('Sensitivity', 'Specificity'), pch=15, col=c('darkorange1', 'cyan4'))

dev.off()


pdf('trainingSet_predictionStats.pdf', width=4, height=5)
barplot(as.matrix(t(trainingData[c(3,2,1),"balanced accuracy"]))
        ,beside=T
        , ylab='Balanced accuracy'
        , ylim=c(0,1)
        , xlab='ATAC peak dynamics'
        ,col = c('darkorchid4','grey', 'goldenrod2')
        ,names.arg = c("Embryo", "No change", "Larval")
)
abline(h=0.5, col='red')

barplot(as.matrix(t(trainingData[c(3,2,1),c("sensitivity", "specificity")]))
        ,beside=T
        , ylab='Balanced accuracy'
        , ylim=c(0,1)
        , xlab='ATAC peak dynamics'
        ,col = c('darkorange1', 'cyan4')
        ,names.arg = c("Embryo", "No change", "Larval")
)
legend('top', legend=c('Sensitivity', 'Specificity'), pch=15, col=c('darkorange1', 'cyan4'))

dev.off()

# I also want to try to combine the training and testing and plot all in one
trainingPlot <- as.matrix(t(trainingData[c(3,2,1),c("sensitivity", "specificity")]))
testingPlot <- as.matrix(t(testingData[c(3,2,1),c("sensitivity", "specificity")]))
combinedPlot <- rbind(testingPlot[1,], trainingPlot[1,], testingPlot[2,], trainingPlot[2,])
row.names(combinedPlot) <- c("Testing sensitivity", "Training sensitivity", "Testing specificity", "Training specificity")

pdf('trainingAndTestingSpecifAndSensiv.pdf', width=6, height=5)
barx <- barplot(combinedPlot
        ,beside=T
        , ylim=c(0,1)
        , xlab='ATAC peak dynamics'
        , ylab='Sensitivity or Specificity'
        ,col = c(rep('darkorange1', 2), rep('cyan4',2))
        , density = rep(c(NA, 30), 2)
        ,names.arg = c("Embryo", "No change", "Larval")
)
legend('top', col=NA
       , c("Sensitivity", "Specificity", "Held out testing set", "Training set")
       , fill=c('darkorange1', 'cyan4', 'black', 'black')
       , pch=15
       , density = c(rep(NA,3), 30)
)

dev.off()



