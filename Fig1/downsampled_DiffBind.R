library(DiffBind)
library(pheatmap)
library(pvclust)
library(sva)

setwd ("/Users/acd13/Desktop/ATAC/Analysis/differential/N2_dev/allButL1/downsampled")
peaks <- read.csv("peaks.bed", header=F, sep="\t")[,1:3]
# Load in Data
all_gDNAMaksedMetaPeaks.gdnaNoL1s <- dba(sampleSheet="downsampled.csv")

countsNoL1s<-dba.count(all_gDNAMaksedMetaPeaks.gdnaNoL1s,peaks=peaks, bParallel=TRUE, score=DBA_SCORE_RPKM_FOLD, bRemoveDuplicates=FALSE,  bScaleControl=TRUE, bCorPlot=TRUE, fragment=1)
# Unfortunately it became apparent that I needed to do some batch removal.  That was done here
# get just the score data
scoreDataNoL1s <-data.frame(countsNoL1s$allvectors[,4:12])
# get the pheontype data
phenoDataNoL1s <- countsNoL1s$samples[,c(2,4,5)]
rownames(phenoDataNoL1s) <- countsNoL1s$samples[,1]
# for combat:
batchNoL1s = phenoDataNoL1s$Replicate
modcombatNoL1s = model.matrix(~1, data=phenoDataNoL1s)
combat_scoreDataNoL1s = ComBat(dat=scoreDataNoL1s, batch=batchNoL1s, mod=modcombatNoL1s, numCovs=NULL, par.prior=TRUE, prior.plots=FALSE)

# Now return adjusted values to diffbind
countsNoL1sBackUpNotComabtted <- countsNoL1s
countsNoL1s$allvectors[,4:12] <- countsNoL1s$vectors[,4:12] <- combat_scoreDataNoL1s

data.count.contrast.withGDNA = dba.contrast(countsNoL1s, minMembers=2, categories=DBA_CONDITION) # we need to tell DiffBind which fall into which groups
data.count.contrast.analyze.withGDNA = dba.analyze(data.count.contrast.withGDNA, bCorPlot=F, 
                                          bParallel=T, bTagwise=F, bFullLibrarySize=T, 
                                          method=DBA_EDGER, bSubControl=T, bReduceObjects=F) # now finally run the main differential analysis


options(max.print=1e8) 
qVal=0.05
# For some reason when looping this wasn't working, so I'm just doing it by hand
i=1
name1<-data.count.contrast.withGDNA$contrasts[[i]]$name1
name2<-data.count.contrast.withGDNA$contrasts[[i]]$name2
Report = dba.report(data.count.contrast.analyze.withGDNA, bCounts=TRUE, bCalledDetail=TRUE, th=qVal, DataType=DBA_DATA_FRAME, contrast=i, method=DBA_EDGER)
sink(file=paste(name1,"_vs_",name2,"_ComBatCorrected_EdgeR_q",qVal,".txt",sep=""))
Report
sink()
i=i+1 #2
name1<-data.count.contrast.withGDNA$contrasts[[i]]$name1
name2<-data.count.contrast.withGDNA$contrasts[[i]]$name2
Report = dba.report(data.count.contrast.analyze.withGDNA, bCounts=TRUE, bCalledDetail=TRUE, th=qVal, DataType=DBA_DATA_FRAME, contrast=i, method=DBA_EDGER)
sink(file=paste(name1,"_vs_",name2,"_ComBatCorrected_EdgeR_q",qVal,".txt",sep=""))
Report
sink()
i=i+2 #4
name1<-data.count.contrast.withGDNA$contrasts[[i]]$name1
name2<-data.count.contrast.withGDNA$contrasts[[i]]$name2
Report = dba.report(data.count.contrast.analyze.withGDNA, bCounts=TRUE, bCalledDetail=TRUE, th=qVal, DataType=DBA_DATA_FRAME, contrast=i, method=DBA_EDGER)
sink(file=paste(name1,"_vs_",name2,"_ComBatCorrected_EdgeR_q",qVal,".txt",sep=""))
Report
sink()
