library(Gviz)
library(GenomicRanges)
library(scales)

Sys.setenv(PATH = "$PATH:/usr/bin:/usr/local/bin:/Users/acd13/Softwares/bedtools2-2.21.0/bin/")


# this is the one line of code that is unique for my plotting.
# all of my enrichment data is log2, but I don't want it to be log2,
# so here I just exponetiate it

properFormattingOfGenomicData <- function(df){
  df$start=as.numeric(as.character(df$start))
  df$end = as.numeric(as.character(df$end))
  if("score" %in% colnames(df)){
    df$score = (2^as.numeric(as.character(df$score)))-1
  }
  return(df)
}


plotAllStagesAndMetaPeaks <- function(locus, windowSize=25, yLower=0, yUpper=NA, plotName=NA, revStrand = FALSE, withChromStates = FALSE, includePreviousEnhancers = FALSE){
  firstSplit <- strsplit(locus, ":", fixed=T)[[1]]
  chrNum <- firstSplit[1]
  
  chrom <- chr <- paste0('chr',chrNum)
  from <- as.numeric(gsub(",","", strsplit(firstSplit[2], "-", fixed=T)[[1]][1]))
  to <- as.numeric(gsub(",","", strsplit(firstSplit[2], "-", fixed=T)[[1]][2]))
  forBedFile <- paste(chrom, from, to, sep="\t")
  plottingRegionBedFile <- "forGviz.bed"
  write.table(forBedFile, file=plottingRegionBedFile, quote=FALSE, row.names=FALSE, col.names=FALSE)
  
  # Add in the peaks
  metaPeaksFile <- "/Users/acd13/Desktop/ATAC/macs2Peaks/usingInserts_75bpShift/N2_dev/withAllReps/allRepsAllButL1MetaPeaks_smartMerged300bpSummitDist.bed"
  metaPeaks <- readInBedFileAsGRange(metaPeaksFile, plottingRegionBedFile)
  if(! is.na(metaPeaks)){
	genome(metaPeaks) = 'ce10'
  	gen <- genome(metaPeaks)
  	metaPeaksTrack <- AnnotationTrack(metaPeaks, name = "ATAC-seq Peaks", col='burlywood4', fill='burlywood4', size=1)
  }
  
  gtrack <- GenomeAxisTrack(fontsize=16)
  gtfFile <- "/Users/acd13/Desktop/geneDefinitions/ce10RefSeqGenesFromIGV.gtf"
  #gtfFile <- "/Users/acd13/Desktop/ce10_RefSeqGenes_10Feb2015.gtf"
  myGeneModels <- makeMyGeneModels(gtfFile, plottingRegionBedFile)
  grtrack <- GeneRegionTrack(myGeneModels, genome = gen, chromosome = chr, name = "Gene Model", size=0.5)
  
  # The commented out code here is if I want to plot the pooled before gDNA normalization results
  baseDir <- "/Users/acd13/Desktop/ATAC/insertSites/N2dev/enrichmentOverGDNA/"
  #EE_atacFile <- paste0(baseDir, "EE_allReps.adjusted.insertSites.",windowSize,"bpWindowMean.log2EnrichOverGDNA.bg.gz")
  #EE_signal <- readInBedGraphFileAsGRange(EE_atacFile, plottingRegionBedFile)

  #L3_atacFile <- paste0(baseDir, "L3_allReps.adjusted.insertSites.",windowSize,"bpWindowMean.log2EnrichOverGDNA.bg.gz")
  #L3_signal <- readInBedGraphFileAsGRange(L3_atacFile, plottingRegionBedFile)

  #YA_atacFile <- paste0(baseDir, "YA_allReps.adjusted.insertSites.",windowSize,"bpWindowMean.log2EnrichOverGDNA.bg.gz")
  #YA_signal <- readInBedGraphFileAsGRange(YA_atacFile, plottingRegionBedFile)

  allReps_signalFile <- paste0("/Users/acd13/Desktop/ATAC/insertSites/N2dev/enrichmentOverGDNA/EE_L3_YA_IndividRepsABC.adjusted.insertSites.singleBP.log2EnrichOverGDNA.meanIn",windowSize, "bpWindows.bed.gz")
  signals <- readInBedGraphFileAsGRangeFromAllReps_combinedStages(allReps_signalFile, plottingRegionBedFile)
  EE_signal <- signals[[1]]
  genome(EE_signal) = 'ce10'
  EE_signalTrack <- DataTrack(EE_signal, name = "Early embryo", fontsize=16, fill.mountain=rep('darkorchid4',2), col='darkorchid4')
  L3_signal <- signals[[2]]
  genome(L3_signal) = 'ce10'
  L3_signalTrack <- DataTrack(L3_signal, name = "Larval Stage 3", fontsize=16)
  YA_signal <- signals[[3]]
  genome(YA_signal) = 'ce10'
  YA_signalTrack <- DataTrack(YA_signal, name = "Young adult", fontsize=16, fill.mountain=rep('darkgreen',2), col='darkgreen')
  
  if(is.na(yUpper)){
    ylims=extendrange(range(c(values(YA_signalTrack), values(EE_signalTrack), values(L3_signalTrack))))
    ylimit <- c(yLower, ylims[2]) 
  }else{
    ylimit <- c(yLower, yUpper)
  }
  
  
  prevEnhTrack = NA
  if(includePreviousEnhancers){
  	# Add in the enhnacers
	prevEnhFile <- "/Users/acd13/Desktop/ATAC/Analysis/knownEnhancers/previouslyDescribedEnhs/approximateLoci.bed"
  	prevEnh <- readInBedFileAsGRange(prevEnhFile, plottingRegionBedFile)
  	if(! is.na(prevEnh)){
		genome(prevEnh) = 'ce10'
  		gen <- genome(prevEnh)
  		prevEnhTrack <- AnnotationTrack(prevEnh, name = "prevEnh", col='red', fill='red', size=1)
  	}
  }
  
  # finally add in the predicted chromatin states
  if(withChromStates){
	  eeChromTrack <- getChromState('EE', plottingRegionBedFile)
  	  l3ChromTrack <- getChromState('L3', plottingRegionBedFile)
  	  yaChromTrack <- getChromState('YA', plottingRegionBedFile)
  	  if(is.na(metaPeaks)){
  	  	if(is.na(prevEnhTrack)){  	  	
	      tracksToPlot <- list(EE_signalTrack, L3_signalTrack, YA_signalTrack, eeChromTrack, l3ChromTrack, yaChromTrack, grtrack, gtrack)
	    }else{
   	      tracksToPlot <- list(EE_signalTrack, L3_signalTrack, YA_signalTrack, prevEnhTrack, eeChromTrack, l3ChromTrack, yaChromTrack, grtrack, gtrack)
   	 	}
	  }else{
  	  	if(is.na(prevEnhTrack)){  	  	
  	      tracksToPlot <- list(EE_signalTrack, L3_signalTrack, YA_signalTrack, metaPeaksTrack, eeChromTrack, l3ChromTrack, yaChromTrack, grtrack, gtrack)
  	  	}else{
  	      tracksToPlot <- list(EE_signalTrack, L3_signalTrack, YA_signalTrack, metaPeaksTrack, prevEnhTrack, eeChromTrack, l3ChromTrack, yaChromTrack, grtrack, gtrack)  	  	
	  	}
	  }
  }else{
  	  if(is.na(metaPeaks)){
	  	if(is.na(prevEnhTrack)){  	  	
	      tracksToPlot <- list(EE_signalTrack, L3_signalTrack, YA_signalTrack, grtrack, gtrack)
	    }else{
   	      tracksToPlot <- list(EE_signalTrack, L3_signalTrack, YA_signalTrack, grtrack, gtrack)
		}
	  }else{
  	  	if(is.na(prevEnhTrack)){  	  	
  	      tracksToPlot <- list(EE_signalTrack, L3_signalTrack, YA_signalTrack, metaPeaksTrack, grtrack, gtrack)
  	  	}else{
  	      tracksToPlot <- list(EE_signalTrack, L3_signalTrack, YA_signalTrack, metaPeaksTrack, prevEnhTrack, grtrack, gtrack)
  	  	}
	  }
  }
  
  if(! is.na(plotName)){pdf(paste0(plotName,'.pdf'), width=4, height=3) }
  plotTracks(tracksToPlot
             , from=from
             , to=to
             , ylim=ylimit
             , add53=TRUE
             ,exponent=3
             ,background.title = NA
             ,col.title="black"
             ,col.axis='black'
             ,reverseStrand=revStrand
  )
  if(! is.na(plotName)){dev.off()}
}


plotEE_v_L3_withDiffPeaks <- function(locus, windowSize=25, comparisonDirection, yLower=0, yUpper=NA, plotName=NA, revStrand = FALSE, withChromStates=FALSE){
  firstSplit <- strsplit(locus, ":", fixed=T)[[1]]
  chrNum <- firstSplit[1]
  
  chrom <- chr <- paste0('chr',chrNum)
  from <- as.numeric(gsub(",","", strsplit(firstSplit[2], "-", fixed=T)[[1]][1]))
  to <- as.numeric(gsub(",","", strsplit(firstSplit[2], "-", fixed=T)[[1]][2]))
  forBedFile <- paste(chrom, from, to, sep="\t")
  plottingRegionBedFile <- "forGviz.bed"
  write.table(forBedFile, file=plottingRegionBedFile, quote=FALSE, row.names=FALSE, col.names=FALSE)
  
  # Add in the peaks
  if(comparisonDirection=="up"){
        peaksFile <- "/Users/acd13/Desktop/ATAC/Analysis/differential/N2_dev/allButL1/differentialReports/combatCorrected/fdr0.05/splitByDirection/EE_vs_L3_ComBatCorrected_EdgeR_q0.05_up.narrowPeak"
        trackName <- "More accessible in embryo"
        trackColor <- "darkorchid4"
  }else if(comparisonDirection=="down"){
        peaksFile <- "/Users/acd13/Desktop/ATAC/Analysis/differential/N2_dev/allButL1/differentialReports/combatCorrected/fdr0.05/splitByDirection/EE_vs_L3_ComBatCorrected_EdgeR_q0.05_down.narrowPeak"
        trackName <- "More accessible in larval"
        trackColor <- "goldenrod2"
  }else{
       warning("unrecognized comparisonDirection. Quitting.")
       return()
  }
  peaks <- readInBedFileAsGRange(peaksFile, plottingRegionBedFile)
  genome(peaks) = 'ce10'
  gen <- genome(peaks)
  peaksTrack <- AnnotationTrack(peaks, name = trackName, col=trackColor, fill=trackColor, size=1)
  
  gtrack <- GenomeAxisTrack(fontsize=16)
  gtfFile <- "/Users/acd13/Desktop/geneDefinitions/ce10RefSeqGenesFromIGV.gtf"
  #gtfFile <- "/Users/acd13/Desktop/ce10_RefSeqGenes_10Feb2015.gtf"
  myGeneModels <- makeMyGeneModels(gtfFile, plottingRegionBedFile)
  grtrack <- GeneRegionTrack(myGeneModels, genome = gen, chromosome = chr, name = "Gene Model", size=0.5)
  
  # I updated this to include all of the replicates
  allReps_signalFile <- paste0("/Users/acd13/Desktop/ATAC/insertSites/N2dev/enrichmentOverGDNA/EE_L3_YA_IndividRepsABC.adjusted.insertSites.singleBP.log2EnrichOverGDNA.meanIn",windowSize, "bpWindows.bed.gz")
  signals <- readInBedGraphFileAsGRangeFromAllReps_combinedStages(allReps_signalFile, plottingRegionBedFile)
  EE_signal <- signals[[1]]
  genome(EE_signal) = 'ce10'
  EE_signalTrack <- DataTrack(EE_signal, name = "Early embryo", fontsize=16, fill.mountain=rep('darkorchid4',2), col='darkorchid4')
  L3_signal <- signals[[2]]
  genome(L3_signal) = 'ce10'
  L3_signalTrack <- DataTrack(L3_signal, name = "Larval Stage 3", fontsize=16)
  
  if(is.na(yUpper)){
    ylims=extendrange(range(c(values(EE_signalTrack), values(L3_signalTrack))))
    ylimit <- c(yLower, ylims[2]) 
  }else{
    ylimit <- c(yLower, yUpper)
  }
  
  # finally add in the predicted chromatin states
  if(withChromStates){
	  eeChromTrack <- getChromState('EE', plottingRegionBedFile)
  	  l3ChromTrack <- getChromState('L3', plottingRegionBedFile)
      tracksToPlot <- list(EE_signalTrack, L3_signalTrack, peaksTrack, eeChromTrack, l3ChromTrack, grtrack, gtrack)
  }else{
      tracksToPlot <- list(EE_signalTrack, L3_signalTrack, peaksTrack, grtrack, gtrack)
  }
  
  if(! is.na(plotName)){pdf(paste0(plotName,'.pdf'), width=4, height=3) }
  plotTracks(tracksToPlot
             , from=from
             , to=to
             , ylim=ylimit
             , add53=TRUE
             ,exponent=3
             ,background.title = NA
             ,col.title="black"
             ,col.axis='black'
             ,reverseStrand=revStrand
  )
  if(! is.na(plotName)){dev.off()}
}

plotL3_v_YA_withDiffPeaks <- function(locus, windowSize=25, comparisonDirection, yLower=0, yUpper=NA, plotName=NA, revStrand = FALSE, withChromStates=FALSE){
  firstSplit <- strsplit(locus, ":", fixed=T)[[1]]
  chrNum <- firstSplit[1]
  
  chrom <- chr <- paste0('chr',chrNum)
  from <- as.numeric(gsub(",","", strsplit(firstSplit[2], "-", fixed=T)[[1]][1]))
  to <- as.numeric(gsub(",","", strsplit(firstSplit[2], "-", fixed=T)[[1]][2]))
  forBedFile <- paste(chrom, from, to, sep="\t")
  plottingRegionBedFile <- "forGviz.bed"
  write.table(forBedFile, file=plottingRegionBedFile, quote=FALSE, row.names=FALSE, col.names=FALSE)
  
  # Add in the peaks
  if(comparisonDirection=="up"){
        peaksFile <- "/Users/acd13/Desktop/ATAC/Analysis/differential/N2_dev/allButL1/differentialReports/combatCorrected/fdr0.05/splitByDirection/L3_vs_YA_ComBatCorrected_EdgeR_q0.05_up.narrowPeak"
        trackName <- "More accessible in larval"
        trackColor <- "goldenrod2"
  }else if(comparisonDirection=="down"){
        peaksFile <- "/Users/acd13/Desktop/ATAC/Analysis/differential/N2_dev/allButL1/differentialReports/combatCorrected/fdr0.05/splitByDirection/L3_vs_YA_ComBatCorrected_EdgeR_q0.05_down.narrowPeak"
        trackName <- "More accessible in adult"
        trackColor <- "darkgreen"
  }else{
       warning("unrecognized comparisonDirection. Quitting.")
       return()
  }
  peaks <- readInBedFileAsGRange(peaksFile, plottingRegionBedFile)
  genome(peaks) = 'ce10'
  gen <- genome(peaks)
  peaksTrack <- AnnotationTrack(peaks, name = trackName, col=trackColor, fill=trackColor, size=1)
  
  gtrack <- GenomeAxisTrack(fontsize=16)
  gtfFile <- "/Users/acd13/Desktop/geneDefinitions/ce10RefSeqGenesFromIGV.gtf"
  #gtfFile <- "/Users/acd13/Desktop/ce10_RefSeqGenes_10Feb2015.gtf"
  myGeneModels <- makeMyGeneModels(gtfFile, plottingRegionBedFile)
  grtrack <- GeneRegionTrack(myGeneModels, genome = gen, chromosome = chr, name = "Gene Model", size=0.5, stacking='dense')
  
  # I updated this to include all of the replicates
  allReps_signalFile <- paste0("/Users/acd13/Desktop/ATAC/insertSites/N2dev/enrichmentOverGDNA/EE_L3_YA_IndividRepsABC.adjusted.insertSites.singleBP.log2EnrichOverGDNA.meanIn",windowSize, "bpWindows.bed.gz")
  signals <- readInBedGraphFileAsGRangeFromAllReps_combinedStages(allReps_signalFile, plottingRegionBedFile)
  YA_signal <- signals[[3]]
  genome(YA_signal) = 'ce10'
  YA_signalTrack <- DataTrack(YA_signal, name = "Young adult", fontsize=16, fill.mountain=rep('darkgreen',2), col='darkgreen')
  L3_signal <- signals[[2]]
  genome(L3_signal) = 'ce10'
  L3_signalTrack <- DataTrack(L3_signal, name = "Larval Stage 3", fontsize=16)
  
  if(is.na(yUpper)){
    ylims=extendrange(range(c(values(YA_signalTrack), values(L3_signalTrack))))
    ylimit <- c(yLower, ylims[2]) 
  }else{
    ylimit <- c(yLower, yUpper)
  }
  
  # finally add in the predicted chromatin states
  if(withChromStates){
	  yaChromTrack <- getChromState('YA', plottingRegionBedFile)
  	  l3ChromTrack <- getChromState('L3', plottingRegionBedFile)
      tracksToPlot <- list(L3_signalTrack, YA_signalTrack, peaksTrack, l3ChromTrack, yaChromTrack, grtrack, gtrack)
  }else{
      tracksToPlot <- list(L3_signalTrack, YA_signalTrack, peaksTrack, grtrack, gtrack)
  }
  
  if(! is.na(plotName)){pdf(paste0(plotName,'.pdf'), width=4, height=3) }
  plotTracks(tracksToPlot
             , from=from
             , to=to
             , ylim=ylimit
             , add53=TRUE
             ,exponent=3
             ,background.title = NA
             ,col.title="black"
             ,col.axis='black'
             ,reverseStrand=revStrand
  )
  if(! is.na(plotName)){dev.off()}
}




readAndFilterToRegion <- function (inputFile, regionFile){
  temp <- strsplit(system(paste("intersectBed -wa -a ", inputFile, " -b ",regionFile, sep=""), intern=TRUE ), "\t")
  return(t(data.frame(sapply(temp, unlist))))
}

readInBedGraphFileAsGRangeFromAllReps <- function(bedFile, plottingRegionBedFile, vectorOfColumnInds, vectorOfIDs=NA){
  allData <- list()
  
  if(all(is.na(vectorOfIDs))){vectorOfIDs <- seq(1, length(vectorOfColumnInds))}
  
  allBedData <- data.frame(readAndFilterToRegion(bedFile, plottingRegionBedFile))
  
  for(i in 1:length(vectorOfColumnInds)){
    bedDataForThisCol <- cbind(allBedData[,1:3], rep(vectorOfIDs[i], nrow(allBedData)), allBedData[,vectorOfColumnInds[i]])
    allData[[i]] <- readInBedFileAsGRangeCore(bedDataForThisCol)
  }
  return(allData)
}


readInBedGraphFileAsGRangeFromAllReps_combinedStages <- function(bedFile, plottingRegionBedFile){
  allData <- list()
  
  vectorOfIDs <- c("EE", "L3", "YA")
  vectorOfColumnInds <- list(EE = c(4:6), L3 = c(7:9), YA = c(10:12))
  
  allBedData <- data.frame(readAndFilterToRegion(bedFile, plottingRegionBedFile))

  for(i in 1:length(vectorOfColumnInds)){
    stageData <- as.matrix(allBedData[,vectorOfColumnInds[[i]]])
    stageMean <- vector(mode='numeric', length=nrow(stageData))
    for (j in 1:nrow(stageData)){stageMean[j] <- mean(as.numeric(stageData[j,]))}
    stageMean[which(stageMean < 0)] <- 0
    bedDataForThisCol <- cbind(allBedData[,1:3], rep(vectorOfIDs[i], nrow(allBedData)), stageMean)
    allData[[i]] <- readInBedFileAsGRangeCore(bedDataForThisCol)
  }
  return(allData)
}

readInBedFileAsGRange <- function(bedFile, plottingRegionBedFile){
  bedData <- data.frame(readAndFilterToRegion(bedFile, plottingRegionBedFile))
  return(readInBedFileAsGRangeCore(bedData))
}

readInBedFileAsGRangeCore <- function(bedData){
  sixColNames <- c('chr','start','end','id','score','strand')
  if(nrow(bedData)==1){
    if(ncol(bedData) > 6){
      names(bedData)[1:6] <-sixColNames
    }else {
      names(bedData) <-sixColNames[1:ncol(bedData)]
    }
  }else if(nrow(bedData) > 1){
    if(ncol(bedData) > 6){
      colnames(bedData)[1:6] <-sixColNames
    }else {
      colnames(bedData) <-sixColNames[1:ncol(bedData)]
    }
  }
  bedData <- properFormattingOfGenomicData(bedData)
  
  if(nrow(bedData)==1){
    if(length(bedData) == 3){
      return(with(bedData, GRanges(chr, IRanges(start, end))))
    }else if(length(bedData) == 4) {
      return(with(bedData, GRanges(chr, IRanges(start, end), id=id)))
    }else if(length(bedData) >= 5) {
      return(with(bedData, GRanges(chr, IRanges(start, end), id=id, score=score)))
      #}else if(length(bedData) >= 6) {
      # return(with(bedData, GRanges(chr, IRanges(start, end), id=id, score=score, strand=strand)))
    }
  }else{
    if(ncol(bedData) ==3){
      return(with(bedData, GRanges(chr, IRanges(start, end))))
    }else if(ncol(bedData) == 4) {
      return(with(bedData, GRanges(chr, IRanges(start, end), id=id)))
    }else if(ncol(bedData) >= 5) {
      return(with(bedData, GRanges(chr, IRanges(start, end), id=id, score=score)))
      # }else if(ncol(bedData) >= 6) {
      #  return(with(bedData, GRanges(chr, IRanges(start, end), id=id, score=score, strand=strand)))
    }
  }
  return(NA)
}

readInBedGraphFileAsGRange <- function(bgFile, plottingRegionBedFile){
  bgData <- data.frame(readAndFilterToRegion(bgFile, plottingRegionBedFile))
  colNames <- c('chr','start','end','score')
  if(nrow(bgData)==1){
    names(bgData)[1:4] <- colNames
  }else if(nrow(bgData) > 1){
    colnames(bgData)[1:4] <- colNames
  }
  bgData <- properFormattingOfGenomicData(bgData)
  
  return(with(bgData, GRanges(chr, IRanges(start, end), score=score)))
}

makeMyGeneModels <- function(gtfFile, plottingRegionBedFile){
  gtfData <- data.frame(readAndFilterToRegion(gtfFile, plottingRegionBedFile))
  ids <- strsplit(as.character(gtfData[,9]), ";")
  geneIds <- gsub("gene_id ","",sapply(ids, '[[',1))
  txIds <- gsub(" transcript_id ","",sapply(ids, '[[',2))
  
  myGeneModels <- cbind(gtfData[,c(1,4,5)], (as.numeric(as.character(gtfData[,5]))-as.numeric(as.character(gtfData[,4]))+1), gtfData[,c(7,3)], geneIds, txIds)
  
  colnames(myGeneModels) <- c("chromosome", "start", "end", "width", "strand", "feature","gene","transcript")
  myGeneModels$chromosome = as.character(myGeneModels$chromosome)
  myGeneModels$feature = as.character(myGeneModels$feature)
  myGeneModels$start = as.numeric(as.character(myGeneModels$start))
  myGeneModels$end = as.numeric(as.character(myGeneModels$end))
  myGeneModels$width = as.numeric(as.character(myGeneModels$width))
  
  return(myGeneModels)
}  


getChromState <- function(stage, plottingRegionBedFile){
    if(stage != 'EE' & stage != 'L3' & stage != 'YA'){warning('Unrecognized stage. Quitting.'); return()}
    
     # the program fails if there is a track without a region, so I'll set up an empty list and only add to it if the region has something in it
    tracksToPlot <- list()

    alLStagesFilePrefix <- "/Users/acd13/Desktop/ChromHMM/ce10_ws220_EE_L3_YA/8marks_EE_L3_YA_1e-3_19states/reordered/beds/"

    TxFile <- paste0(alLStagesFilePrefix, stage , "/all", stage, "TxMerged.bed")
    Tx <- readInBedFileAsGRange(TxFile, plottingRegionBedFile)
    if(! is.na(Tx) ){ tracksToPlot[[ length(tracksToPlot) +1 ]] <- AnnotationTrack(Tx, name = "Transcribed", col='palegreen3', fill='palegreen3') }
	
	TssFile <- paste0(alLStagesFilePrefix, stage , "/all", stage, "TSSMerged.bed")
	Tss <- readInBedFileAsGRange(TssFile, plottingRegionBedFile)
	if(! is.na(Tss) ){ tracksToPlot[[ length(tracksToPlot) +1 ]] <- AnnotationTrack(Tss, name = "Promoter", col='tomato3', fill='tomato3') }
	
	ActEnhFile <- paste0(alLStagesFilePrefix, stage , "/", stage, "ActiveEnhancers_notRepEnh13.bed")
	ActEnh <- readInBedFileAsGRange(ActEnhFile, plottingRegionBedFile)
	if(! is.na(ActEnh) ){ tracksToPlot[[ length(tracksToPlot) +1 ]] <- AnnotationTrack(ActEnh, name = "Active Enhancer", col='sienna1', fill='sienna1') }
		
	RepEnhFile <- paste0(alLStagesFilePrefix, stage , "/", stage, "_RepEnh.bed")
	RepEnh <- readInBedFileAsGRange(RepEnhFile, plottingRegionBedFile)
	if(! is.na(RepEnh) ){ tracksToPlot[[ length(tracksToPlot) +1 ]] <- AnnotationTrack(RepEnh, name = "Rep'd Enhnacer", col='darkkhaki', fill='darkkhaki') }

	RepressedRegionsFile <- paste0(alLStagesFilePrefix, stage , "/all", stage, "RepMerged_notRepEnh13.bed")
	RepressedRegions <- readInBedFileAsGRange(RepressedRegionsFile, plottingRegionBedFile)
	if(! is.na(RepressedRegions) ){ tracksToPlot[[ length(tracksToPlot) +1 ]] <- AnnotationTrack(RepressedRegions, name = "Repressed", col='gray65', fill='gray65') }

	HcFile <- paste0(alLStagesFilePrefix, stage , "/all", stage, "HetChromMerged.bed")
	Hc <- readInBedFileAsGRange(HcFile, plottingRegionBedFile)
	if(! is.na(Hc) ){ tracksToPlot[[ length(tracksToPlot) +1 ]] <- AnnotationTrack(Hc, name = "Heterochromatin", col='slateblue3', fill='slateblue3') }

	OthersFile <- paste0(alLStagesFilePrefix, stage , "/", stage, "NotInActiveHetRep.bed")
	OthersRegions <- readInBedFileAsGRange(OthersFile, plottingRegionBedFile)
	if(! is.na(OthersRegions) ){ tracksToPlot[[ length(tracksToPlot) +1 ]] <- AnnotationTrack(OthersRegions, name = "Low/Other", col=NA, fill='gray88') }
  
    
	return( OverlayTrack(tracksToPlot ,background.title = NA) )
}





#=================
# set up a scheme
#===============
scheme <- getScheme()
#scheme$GeneRegionTrack$collapseTranscripts='meta'
scheme$GeneRegionTrack$fill <- 'black'
scheme$GeneRegionTrack$col <- NULL
scheme$GeneRegionTrack$transcriptAnnotation <- NULL
scheme$GeneRegionTrack$shape <- c('box', 'arrow')
scheme$GeneRegionTrack$stacking <- 'dense'
scheme$GeneRegionTrack$col.axis <- 'black'
scheme$GeneRegionTrack$col.title <- "black"

scheme$DataTrack$fill.mountain=rep('goldenrod2',2)
scheme$DataTrack$col='goldenrod2'
scheme$DataTrack$type='polygon'
scheme$DataTrack$col.axis <- 'black'
scheme$DataTrack$col.title <- "black"

scheme$AnnotationTrack$col.title <- "black"
scheme$AnnotationTrack$col.axis <- 'black'
scheme$AnnotationTrack$size <- 0.8
scheme$AnnotationTrack$fontsize <- 12

addScheme(scheme, "myScheme")
options(Gviz.scheme = "myScheme")

