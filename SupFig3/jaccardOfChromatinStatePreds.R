# originally written 13Jul2015 by Aaron Daugherty - Brunet Lab, Stanford University
# this program basically calclulates all pairwise jaccards for a set of genomic locations
# in this case I'm comparing between different chomatin state model predictions

# set up
Sys.setenv(PATH = "$PATH:/usr/bin:/usr/local/bin:/Users/acd13/Softwares/bedtools2-2.21.0/bin/")
library(pheatmap)
library(pvclust)
setwd('/Users/acd13/Desktop/ATAC/Analysis/checkingChromHMMs')

genJacMat <- function(fileListOne, fileListTwo, namesOne, namesTwo){
  jaccards.nucleotides <- matrix(0, nrow=length(fileListOne),ncol=length(fileListTwo))
  rownames(jaccards.nucleotides) <- namesOne
  colnames(jaccards.nucleotides) <- namesTwo
  # run all of the pairwise comparisons by using system calls to Bedtools
  for (i in 1:length(fileListOne)){
    for (j in 1:length(fileListTwo)){
      # bedtools does it counting the number of bp overlap
      rawResult <- system(paste("bedtools jaccard -a ", fileListOne[i], " -b ", fileListTwo[j],sep=""), intern=TRUE)
      jaccards.nucleotides[i,j] <- as.numeric(strsplit(rawResult[2], "\t")[[1]][3])
    }
  }
  return(jaccards.nucleotides)
}


stage='EE'
hiHMMStateDir=paste0(stage, "_hiHMMStates")
# get the lists of my files
fileOfMyPreds <- read.table(paste(stage, "ChromHMMStatesToCompareTo.txt", sep=""), header=F)[,1]
myMergedPredFiles <- fileOfMyPreds[c(1:6, 13)]
myMergedPredNames <- strsplit(as.character(myMergedPredFiles), "/", fixed=T)
myMergedPredNames <- gsub("_", "", gsub(stage, "", gsub(".bed", "", sapply(myMergedPredNames, "[[", length(myMergedPredNames[[1]])))))

myIndividualPredFiles <- fileOfMyPreds[-c(1:6, 13)]
myIndividualPredNames <- strsplit(as.character(myIndividualPredFiles), "/", fixed=T)
myIndividualPredNames <- gsub("_", "", gsub(stage, "", gsub(".bed", "", sapply(myIndividualPredNames, "[[", length(myIndividualPredNames[[1]])))))

# get the list of their files
theirIndividPredFiles <- list.files(hiHMMStateDir, pattern="*.bed$")
theirIndividPredNames <- gsub(".bed", "", theirIndividPredFiles)
# the list files doesn't include the subdirectory, so we add that here
for(i in 1:length(theirIndividPredFiles)){
  theirIndividPredFiles[i] <- paste0(hiHMMStateDir,"/",theirIndividPredFiles[i])
}

theirMergedPredFiles <- list.files(paste0(hiHMMStateDir, "/merged"), pattern="*.bed$")
theirMergedPredNames <- gsub("merged", "", gsub(".bed", "", theirMergedPredFiles))
# the list files doesn't include the subdirectory, so we add that here
for(i in 1:length(theirMergedPredFiles)){
    theirMergedPredFiles[i] <- paste0(hiHMMStateDir,"/merged/",theirMergedPredFiles[i])
}
# they had some single states that I want to include here
# the promoter state
theirMergedPredFiles[length(theirMergedPredFiles) + 1] <- theirIndividPredFiles[1]
theirMergedPredNames[length(theirMergedPredNames) + 1] <- theirIndividPredNames[1]
# next the unmapped state
theirMergedPredFiles[length(theirMergedPredFiles) + 1] <- theirIndividPredFiles[9]
theirMergedPredNames[length(theirMergedPredNames) + 1] <- theirIndividPredNames[9]
# and finally a genic H4K20me1 state
theirMergedPredFiles[length(theirMergedPredFiles) + 1] <- theirIndividPredFiles[14]
theirMergedPredNames[length(theirMergedPredNames) + 1] <- theirIndividPredNames[14]


# I want to do all of the possible comparisons: both my individ and my merged vs both their individ and their merged

myMergedVsAllTheirMerged <- genJacMat(theirMergedPredFiles, myMergedPredFiles, theirMergedPredNames, myMergedPredNames)
myMergedVsTheirCoreMerged <- genJacMat(theirMergedPredFiles[-c(1,2,6, 10, 11)], myMergedPredFiles, theirMergedPredNames[-c(1,2,6, 10, 11)], myMergedPredNames)
#
myMergedVsTheirIndivid <- genJacMat(theirIndividPredFiles, myMergedPredFiles, theirIndividPredNames, myMergedPredNames)
#
myIndividVsTheirMerged <- genJacMat(theirMergedPredFiles, myIndividualPredFiles, theirMergedPredNames, myIndividualPredNames)
myIndividVsTheirCoreMerged <- genJacMat(theirMergedPredFiles[-c(1,2,6, 10, 11)], myIndividualPredFiles, theirMergedPredNames[-c(1,2,6, 10, 11)], myIndividualPredNames)
#
myIndividVsTheirIndivid <- genJacMat(theirIndividPredFiles, myIndividualPredFiles, theirIndividPredNames, myIndividualPredNames)

myMergedVsTheirCoreMerged_forPlotting <- genJacMat(theirMergedPredFiles[c(3,9,4,8,5,7)], myMergedPredFiles[c(6,5,2,7,3,4,1)], theirMergedPredNames[c(3,9,4,8,5,7)], myMergedPredNames[c(6,5,2,7,3,4,1)])


pdf(paste0(stage, "_pheatmaps_chromHMMVsHiHMM.pdf"), width=5, height=5)
  pheatmap(myMergedVsTheirCoreMerged_forPlotting, cluster_rows=F, cluster_cols=F, main=paste0(stage, "_mergedStatesCondensed"))
  pheatmap(myMergedVsAllTheirMerged, main=paste0(stage, "_mergedStates"))
  pheatmap(myMergedVsTheirIndivid, main=paste0(stage, "_myMergedVsTheirIndivid"))
  pheatmap(myIndividVsTheirMerged, main=paste0(stage, "_myIndividVsTheirMerged"))
  pheatmap(myIndividVsTheirIndivid, cluster_rows=F, cluster_cols=F, main=paste0(stage, "_myIndividVsTheirIndivid"))
dev.off()

