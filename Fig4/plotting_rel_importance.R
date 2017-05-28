# show results
  pdf(paste(baseFilehandle, "_gbmCourseTuning.pdf", sep=""), width=5, height=5)
    plot(gbms)
    gbmmodWithTraining <- gbms$finalModel
    importance<- summary(gbmmodWithTraining, n.trees=length(gbmmodWithTraining$trees))
    barx <- barplot(importance[,2], 
                ylab='Relative importance for final regression model', 
                xaxt='n'
    )
  dev.off()