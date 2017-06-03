library(pheatmap)

plotWilcoxResults <- function(listToUse, fileNamesList, pdfName=NA){
  if(! is.na(pdfName)){pdf(pdfName, width=length(listToUse))}
  allNegLog10Qvals <- matrix(0, ncol=length(listToUse), nrow=length(listToUse))
  rownames(allNegLog10Qvals) <- colnames(allNegLog10Qvals) <- fileNamesList
  boxplot(listToUse, notch=T, outline=F,names=fileNamesList, xlab="ATAC-seq fragment size")
  
  for(i in 1:length(listToUse)){
    for(j in 1:length(listToUse)){
      allWilcoxTest <- wilcox.test(listToUse[[i]], listToUse[[j]], alternative='greater') # I'm specifically interested in if perPeakList[[i]] has greater insert size distributions
      allNegLog10Qvals[i,j] <- -log10(as.numeric(gsub("<", "", p.adjust(format.pval(allWilcoxTest$p.value,6,1e-323), method='BH', n=length(listToUse)^2))))
    }
  }
  pheatmap(allNegLog10Qvals)
  if(! is.na(pdfName)){dev.off()}
  return(allNegLog10Qvals)
}

getMedianInsertList <- function(dir){
  allFiles <- list.files(path=dir, pattern="*.txt$")
  
  if( length(allFiles) == 0 ){
    warning(paste("No files matching *.txt found in", dir))
    quit("no", 13)
  }
  
  allData <- list()
  fileNamesList <- vector(mode='character', length=length(allFiles))
  for (i in 1:length(allFiles) ) {
    temp <- unique(read.csv(paste0(dir,"/", allFiles[i]), header=F, sep="\t"))
    allData[[i]] <- as.numeric(as.character(temp[,2]))
    names(allData[[i]]) <- as.character(temp[,1])
    fileNamesList[i] <- gsub(".txt", "",  allFiles[i])
  }
  names(allData) <- fileNamesList
  return(allData)
}

getPerMillionMappedData <- function(dir, totalReadCount){
  
  allFiles <- list.files(path=dir, pattern="*.narrowPeak$") # these aren't actually narrowPeak, but that's what they got called in the processing
  
  if( length(allFiles) == 0 ){
    warning(paste("No .narrowPeak files found in", dir))
    quit("no", 13)
  }
  
  countsData <- list()
  fileNamesList <- vector(mode='character', length=length(allFiles))
  for (i in 1:length(allFiles) ) {
    temp <- unique(read.csv(paste0(dir,"/", allFiles[i]), header=F, sep="\t"))
    countsData[[i]] <- unique(data.frame(peakName = as.character(temp[,1]), pmm = as.numeric(as.character(temp[,2]))/totalReadCount*1e6))
    fileNamesList[i] <- gsub(".narrowPeak", "",  allFiles[i])
  }
  names(countsData) <- fileNamesList
  
  return(countsData)
}

getEnrichmentData <- function(stageCountsDir, stageTotalCounts, controlPMMData){

   stagePMM <- getPerMillionMappedData(stageCountsDir, stageTotalCounts)
   enrichData <- list()

   for(i in 1:length(stagePMM)){

     # this is overly cautious, but this doesn't assume that the control data was saved in the same order
     controlDataIndex <- grep(names(stagePMM)[i], names(controlPMMData))

     # first we merge the 2 datasets, that has the advantage of ensuring they both cover the same points
     originalLength <- nrow(unique(stagePMM[[i]]))
     withControl <- unique(merge(controlPMMData[[ controlDataIndex ]], stagePMM[[i]], by = 1))
     
     if(nrow(withControl) != originalLength){
       warning(paste0("Differing numnber of peaks before and after merge. Before: ", originalLength, ". And after: ", nrow(withControl)))
       quit('no', 13)
     }
     
     enrichData[[i]] <- log2( (withControl[,3] + toAvoidZero) / (withControl[,2] + toAvoidZero )) 
     names(enrichData[[i]]) <- withControl[,1]
  }
  
  names(enrichData) <- names(stagePMM)
  return(enrichData)
}                                

namedVectorToDataFrame <- function(namedVector){
  return(data.frame(Names = names(namedVector), values = namedVector))
}

compareNamedVectors <- function(namedVectorOne, namedVectorTwo){
  firstDF <- namedVectorToDataFrame(namedVectorOne)
  secondDF <- namedVectorToDataFrame(namedVectorTwo)
  
  combinedDF <- unique(merge(firstDF, secondDF, by='Names'))
  
  if(nrow(combinedDF) != nrow(unique(firstDF))){
    warning(paste0("In combining DFs, differing number of rows found. Started with: ", nrow(unique(firstDF)),". Ended with: ",  nrow(combinedDF)))
    quit('no', 13)
  }
  
  toReturn <- combinedDF[,2] - combinedDF[,3]
  names(toReturn) <- combinedDF[,1]
  return(toReturn)
}
  