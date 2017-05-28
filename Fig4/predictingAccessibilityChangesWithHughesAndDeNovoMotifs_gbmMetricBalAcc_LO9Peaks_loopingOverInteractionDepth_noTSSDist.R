# Updates to this version is that we now use de novo mapped motifs, and do not include the distance to the TSS b/c it wasn't helping anyways

library(tools)
library(caret)
library(gbm)
library(doMC)
registerDoMC(cores = 6)

# define functions
# this function reads in the file and parses through the way I encoded the occurence count
getOccurencesIntoDf <- function(occurencesFile, featureDf, regionDescriptor=""){
          motifOccurence <- read.table(occurencesFile, header=F, sep="\t")
          # the format of this file is:
          # motif name <tab> csv of peaks where this motif was seen.
          # each peak is in the following format:
          # chr_start_stop_numberOfOccurencesOfMotif
          
          # I'll just update the cells where there are occurences, and put in the 0s later
          for(i in 1:nrow(motifOccurence)){ # for each motif
              # pull out the peaks list make each peak its own
              allPeaks <- strsplit(as.character(motifOccurence[i,2]), ",", fixed=T)[[1]]
              # the occurence number ended up being encoded with the name, so we have to pull it out, and then put the name back together in a DF
              splitApart <- strsplit(as.character(allPeaks), "_", fixed=T)
              occurences <- sapply(splitApart, "[[", 4)
              nameParts <- cbind(sapply(splitApart, "[[", 1), sapply(splitApart, "[[", 2), sapply(splitApart, "[[", 3))
              names <- apply(nameParts, 1, function(x){paste0(x,collapse="_")})
              
              # take advantage of merge to make this go much faster
              forMerging <- data.frame(name=names, counts=occurences)
              featureDf <- unique(merge(featureDf, forMerging, by=1, all.x=T)) # the all.x keeps all the peaks in the data frame
              colnames(featureDf)[i+1] <- paste(regionDescriptor, as.character(motifOccurence[i,1]), sep="_") # name the column after the motif
          }
          return(featureDf)
}


getMultiClassConfusionMatrixBalancedAccuracy <- function (data, lev = NULL, model = NULL) {
    
    if (!all(levels(data[, "pred"]) == levels(data[, "obs"]))){ 
        stop("levels of observed and predicted data do not match")
    }
    
    testConfMat <- try(confusionMatrix(data$pred, data$obs), silent = TRUE)
    if (class(testConfMat)[1] == "try-error") {
        return(NA)
    }else{
        out <- mean(testConfMat$byClass[,"Balanced Accuracy"])
        names(out) <- "balancedAcc"
        return(out)
    }
}


#=====================================================================================================================
# read in raw data and set up
#=====================================================================================================================

setwd("/Users/acd13/Desktop/ATAC/Analysis/predictingAccessChangeWithMotifs")

peakMotifOccurences <- "./motifCounts/hughesMotifAndDeNovoCountsInMetaAtacPeaks_LO9.motifOccurenceList.txt"

interactionDepthRange = seq(2,11,3)

peakLabelsFile <- "/Users/acd13/Desktop/ATAC/Analysis/predictingAccessChangeWithMotifs/metaAtacPeaks_EE_L3_fdr0.05Labels.txt"

prefix <- paste("usingGBM/predicting", file_path_sans_ext(basename(peakLabelsFile)), "withLO9HughesAndDeNovoMotifs_peaks_noTssDist", sep="_")
dir.create(prefix)

#=====================================================================================================================
# generate feature matrix
#=====================================================================================================================

# we can generate the feature matrix regardless of weighting
outputFH <- paste0(prefix,"/predicting_", file_path_sans_ext(basename(peakLabelsFile)), "_withLO9HughesAndDeNovoMotifs_peaks")
ftMatRds <- paste(outputFH, "featMat.rds", sep="_")

if(file.exists(ftMatRds)){
  
  features <- readRDS(ftMatRds)
  print(paste("Using saved features matrix", ftMatRds))
  
  colNamesOfFeatMat <- read.table(paste(outputFH, "colnamesFM.tsv", sep="_"), sep="\t", header=T)[,1]
  colnames(features)[1:(ncol(features)-1)] <- seq(1,(ncol(features)-1))
  indsToTest <- readRDS(paste(outputFH,"testingInds.rds",sep="_"))
  testSet <- features[indsToTest,]
  fullTrainingSet <- features[-indsToTest,]
  
}else{
  
  peakLabels <- read.table(peakLabelsFile, header=F, sep="\t")
  # the format of this file is:
  # peak name <tab> label (1,0,-1)
  # where peak name is chr_start_stop
  # and the peak label has to do with wether it is differentially accessible and if so in which direction
  
  # i want to turn all of this into a matrix where each row is a peak, and each column is a motif
  # we'll add the motifs below, but first initializ the DF with the peaks
  featureDf <- data.frame(peaks=peakLabels[,1])
  
  # now call a function that reads in the occurences and processes them as necessary
  # we'll provide the same starting df, and then merge all three afterward
  
  
  peakOccurencesDf <- getOccurencesIntoDf(peakMotifOccurences, featureDf, "peaks")
  #
  # Now combine the three
  featureDf <- peakOccurencesDf
  # add in the peak labels  
  featureDf <- unique(merge(featureDf, peakLabels, by=1))
  
  # the data frame and its stupid factors were causing issues, so let's avoid that
  # the important thing being that they stay in the same order, allowing us to bind them back together later
  featureMatrix <- as.matrix(featureDf[,-c(1, ncol(featureDf))])
  # now replace the NAs with 0
  featureMatrix[is.na(featureMatrix)] <- 0
  # and put things back together
  features <- as.data.frame(cbind(featureMatrix, featureDf[,ncol(featureDf)]))
  rownames(features) <- featureDf[,1]
  features[,ncol(features)] <- as.factor(features[,ncol(features)])
  colnames(features)[ncol(features)] <- "peakLabels"
  
  # I want the features to be numeric
  for(k in 1:(ncol(features)-1)) { features[,k] <- as.numeric(as.character(features[,k])) }
  
  saveRDS(features, file=ftMatRds)
  print("Feature matrix generated and saved")
  rm(featureDf, featureMatrix, peakOccurencesDf, peakLabels, allOccurences)
  
  
  # I'll use a subset of the peaks to optimize parameters, then build a model using all the training set, and then finally test on the remaining peaks
  # subsample the peaks
  
  # I was getting some weird things with the column names, so instead of trying to keep the informative names, I'll just have a key that I'll decode later
  colNamesOfFeatMat <- colnames(features)
  write.table(colNamesOfFeatMat, file=paste(outputFH, "colnamesFM.tsv", sep="_"), quote=F, sep="\t", row.names=F)
  colnames(features)[1:(ncol(features)-1)] <- seq(1,(ncol(features)-1))
  
  # now back to work
  numberToTest <-floor(nrow(features)*(1/3))
  indsToTest <- sample(1:nrow(features), numberToTest)
  testSet <- features[indsToTest,]
  
  fullTrainingSet <- features[-indsToTest,]
  saveRDS(indsToTest, file=paste(outputFH,"testingInds.rds",sep="_"))
}

outputFHOrig <- outputFH
#=====================================================================================================================
# Generate the GBM
#=====================================================================================================================
for(interactionDepth in interactionDepthRange){
    print(paste0("Starting interactionDepth ", interactionDepth))
    outputFH <- paste(outputFHOrig, "interactionDepth", interactionDepth, "metricBalAcc", sep="_")

    if(file.exists(paste(outputFH, "featureImportance.tsv", sep="_"))){next}
    
    gbmOutFile <- paste(outputFH, "caretdGBM.rds", sep="_")
  
    if(file.exists(gbmOutFile)){
        myOptdGbm <- readRDS(gbmOutFile)
        print(paste("Using saved GBM", gbmOutFile))
    }else{
         
        # set up the option controls for caret
        my.ctrl.opt <-trainControl(method = "cv", number = 10, classProbs=TRUE, allowParallel=TRUE, verbose=T, summaryFunction=getMultiClassConfusionMatrixBalancedAccuracy)
        courseGbmGrid <- expand.grid(.interaction.depth = interactionDepth, .n.trees = seq(6000,14000, 4000), .shrinkage = 0.001)
        
        print("Print beginning to train GBM")
        myOptdGbm <- train(x = fullTrainingSet[,-ncol(fullTrainingSet)], y = fullTrainingSet[,ncol(fullTrainingSet)]
                      , method = "gbm", trControl = my.ctrl.opt
                      , tuneGrid = courseGbmGrid, metric="balancedAcc", maximize=TRUE, verbose=T)
        saveRDS(myOptdGbm, file = gbmOutFile)
    
        print("GBM finished and saved. Next printing out performance stats")
    }
    #=====================================================================================================================
    # Make predictions on the test set and calculate metrics
    #=====================================================================================================================
    
    testingPr <- predict(myOptdGbm$finalModel, type='response', newdata=testSet[,-ncol(testSet)], n.trees= myOptdGbm$finalModel$n.trees)[,,1]
    assignmentPreds <- apply(testingPr, 1, function(x){ colnames(testingPr)[which.max(x)]})
    forConfMat <- as.factor(c(seq(-1,1),assignmentPreds)) # when I was testing this with very few trees it sometimes didn't predict anything for a class, so I include these, just for safety
    testConfMat <- confusionMatrix(forConfMat[-c(1:3)], testSet[,ncol(testSet)])
    for(ind in 1:length(testConfMat)){write.table(testConfMat[[ind]], file=paste(outputFH, "testSetConfusionMatrixAndStats.txt", sep="_"), append=TRUE, quote=F, sep="\t", col.names=F)}
    
    pdf(paste0(outputFH, "_gbmClassificationPlots.pdf"))
        plot(myOptdGbm)
        importance <- summary(myOptdGbm$finalModel, n.trees=myOptdGbm$finalModel$n.trees)
    dev.off()
    write.table(importance, file=paste(outputFH, "featureImportance.tsv", sep="_"), quote=F, sep="\t")
    
    # I also need to pull out how the model did on the training set
    trainingAssignments <- apply(myOptdGbm$finalModel$fit, 1, function(x){ colnames(myOptdGbm$finalModel$fit)[which.max(x)]})
    forConfMat <- as.factor(c(seq(-1,1),trainingAssignments)) # when I was testing this with very few trees it sometimes didn't predict anything for a class, so I include these, just for safety
    trainConfMat <- confusionMatrix(forConfMat[-c(1:3)], fullTrainingSet[,ncol(fullTrainingSet)])
    for(ind in 1:length(trainConfMat)){write.table(trainConfMat[[ind]], file=paste(outputFH, "trainSetConfusionMatrixAndStats.txt", sep="_"), append=TRUE, quote=F, sep="\t", col.names=F)}
}


