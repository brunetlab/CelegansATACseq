library(pheatmap)
library(RColorBrewer)
library(R.utils)
getStageFileNames <- function(stage){
  if(stage=='L3'){
    stageChromStateNames <- c(paste0("all", stage, "TxMerged")
                              , paste0("all", stage, "TSSMerged")
                              , paste0(stage, "_allActiveEnhancers_notRepEnh13")
                              , paste0(stage, "_repressedEnh")
                              , paste0("all", stage, "RepMerged_notRepEnh13")
                              , paste0("all", stage, "HetChromMerged")
                              , paste0(stage, "NotInActiveHetRep")
    )
  }else if (stage=='EE'){
    stageChromStateNames <- c(paste0("all", stage, "TxMerged")
                              , paste0("all", stage, "TSSMerged")
                              , paste0(stage, "ActiveEnhancers_notRepEnh13")
                              , paste0(stage, "_RepEnh")
                              , paste0("all", stage, "RepMerged_notRepEnh13")
                              , paste0("all", stage, "HetChromMerged")
                              , paste0(stage, "NotInActiveHetRep")
    )
  }else if (stage=='YA'){
    stageChromStateNames <- c(paste0("all", stage, "TxMerged")
                              , paste0("all", stage, "TSSMerged")
                              , paste0(stage, "ActiveEnhancers_notRepEnh13")
                              , paste0(stage, "_RepEnh")
                              , paste0("all", stage, "RepMerged_notRepEnh13")
                              , paste0("all", stage, "HetChromMerged")
                              , paste0(stage, "NotInActiveHetRep")
    )
  }
  return(stageChromStateNames)
}



stage1="L3"
stage2="YA"

mainDir <- paste("/Users/acd13/Desktop/ATAC/Analysis/changingChromStates/annotBedFile", stage1, "vs", stage2, sep="_")
setwd(mainDir)

baseName <- paste(stage1, "v", stage2, "fdr0.05", sep="_")
directionalityNames <- c("upPeaks", "downPeaks", "noChangePeaks")
firstStageChromStateNames <- getStageFileNames(stage1)
secondStageChromStateNames <- getStageFileNames(stage2)

if(file.exists(paste(baseName,"chromHMMStateChanges.rds", sep="_"))){
  allMatrices <- readRDS(paste(baseName,"chromHMMStateChanges.rds", sep="_"))
}else{
  
  allMatrices <- list()
  
  # none of this is so complex, so we'll just nest for loops for ease
  for (k in 1:length(directionalityNames)){
    direction <- directionalityNames[k]
    
    # create the matrix that is going to hold the results
    # the extra space is to note the number of peaks that were lost or unnanotated (which would happen if they were <50 in any one state)
    counts <- matrix(0, nrow=(length(firstStageChromStateNames)+1), ncol=(length(secondStageChromStateNames)+1))
    colnames(counts) <- c(secondStageChromStateNames, "unannotated")
    rownames(counts) <- c(firstStageChromStateNames, "unannotated")
    
    prefix = paste(baseName, direction, sep="_")
    # so we know how many peaks we started with for this direction
    directionPeakCount = countLines(paste(prefix,".bed", sep=""))
    
    firstStageChromCounts <- vector(mode='numeric',length=length(firstStageChromStateNames))
    
    for(i in 1:length(firstStageChromStateNames)){
      firstStageChromPrefix = paste(prefix, firstStageChromStateNames[i], sep="_")
      # so we know how many peaks we started with for this particular set of peaks
      firstStageChromCounts[i] = countLines(paste(firstStageChromPrefix,".bed", sep=""))
      
      for(j in 1:length(secondStageChromStateNames)){
        counts[i,j] <- countLines(paste(firstStageChromPrefix,"_", secondStageChromStateNames[j], ".bed", sep=""))
      }
      # figure out how many were unnanotated (see above)
      counts[i,ncol(counts)] <- firstStageChromCounts[i] - sum(as.numeric(counts[i,1:(ncol(counts)-1)]))
    }
    counts[nrow(counts),1] <- directionPeakCount - sum(firstStageChromCounts)
    
    allMatrices[[k]] <- counts
  }
  
  names(allMatrices) <- directionalityNames
  
  saveRDS(allMatrices, file=paste(baseName,"chromHMMStateChanges.rds", sep="_"))
}

# with the data read in, now let's normalize things to a portion
# this is done on a row-by-row basis to allow for comparison
# this is a little inefficient in that I could have included this above, but this is just easier
rowNormdMats <- vector(mode='list', length=length(directionalityNames)) 
for (k in 1:length(directionalityNames)){
  # in this case I don't care about the peaks that weren't annotated in the first stage
  # but I will keep the ones that were lost in the transition, mostly out of curiosity
  rowNormdMats[[k]] <- matrix(0, nrow=length(firstStageChromStateNames), ncol=length(secondStageChromStateNames)+1) 
  for(i in 1:length(firstStageChromStateNames)){
    rowNormdMats[[k]][i,] <- allMatrices[[k]][i,]/sum(allMatrices[[k]][i,])
  }
}

# now we can compare each set of direction change with the no change peaks (as those are sort of a background to see what happens with other peaks)
# there are some 0s, so to avoid inf values, I'll add a small value
# rather than arbitraily choosing one, I'll use half the smallest-non-0 in this data
toAddToAvoidZero <- min(unlist(rowNormdMats)[which(unlist(rowNormdMats) > 0)])/2
upPeakEnrichVsNoChange <- log2( (rowNormdMats[[1]] + toAddToAvoidZero ) / (rowNormdMats[[3]] + toAddToAvoidZero) )
downPeakEnrichVsNoChange <- log2( (rowNormdMats[[2]] + toAddToAvoidZero ) / (rowNormdMats[[3]] + toAddToAvoidZero) )

colnames(upPeakEnrichVsNoChange) <- colnames(downPeakEnrichVsNoChange) <- c(paste(stage2, "Transcribed")
                                                                            , paste(stage2, " TSS/Promoter")
                                                                            , paste(stage2, " Active Enhancer")
                                                                            , paste(stage2, " Repressed Enhancer")
                                                                            , paste(stage2, " H3K27me3 Repressed")
                                                                            , paste(stage2, " Heterochromatin")
                                                                            , paste(stage2, " Low")
                                                                            , "unannotated"
                                                                            )
rownames(upPeakEnrichVsNoChange) <- rownames(downPeakEnrichVsNoChange) <- c(paste(stage1, "Transcribed")
                                                                               , paste(stage1, " TSS/Promoter")
                                                                               , paste(stage1, " Active Enhancer")
                                                                               , paste(stage1, " Repressed Enhancer")
                                                                               , paste(stage1, " H3K27me3 Repressed")
                                                                               , paste(stage1, " Heterochromatin")
                                                                               , paste(stage1, " Low")
                                                                            )

pdf(paste("../chromatinStateChangeEnrichments_in_", baseName, "_RedBlueHeatmap.pdf", sep=""), width=6, height=6)
  pheatmap(upPeakEnrichVsNoChange[,-ncol(upPeakEnrichVsNoChange)]
           , cluster_rows=F, cluster_cols=F
           , main='up'
           ,  color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100))
  pheatmap(downPeakEnrichVsNoChange[,-ncol(downPeakEnrichVsNoChange)]
         , cluster_rows=F, cluster_cols=F
         , main='down'
         ,  color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100))
dev.off()

