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

combatCorNoL1s <- cor(combat_scoreDataNoL1s)
combatCorSpearmanNoL1s <- cor(combat_scoreDataNoL1s, method='spearman')
pvlusted.combatCorNoL1s <- pvclust(combatCorNoL1s, method.hclust='complete')
pvlusted.combatCorSpearmanNoL1s <- pvclust(combatCorSpearmanNoL1s, method.hclust='complete')

pdf("noL1s_DBA_SCORE_RPKM_FOLD_combated.pdf")

plot(pvlusted.combatCorNoL1s)
pheatmap(combatCorNoL1s)
plot(pvlusted.combatCorSpearmanNoL1s)
pheatmap(combatCorSpearmanNoL1s)

dev.off()

# Now return adjusted values to diffbind
countsNoL1sBackUpNotComabtted <- countsNoL1s
countsNoL1s$allvectors[,4:12] <- countsNoL1s$vectors[,4:12] <- combat_scoreDataNoL1s


inPeakCors.usingGdna<-dba.plotHeatmap(countsNoL1s,minval=0,maxval=1,colScheme="Blues",ColAttributes=NULL)


pdf("CountsInPeaks_Heatmap_subControl_ComBatCorrected.pdf") 
dba.plotHeatmap(countsNoL1s,minval=0,maxval=1,colScheme="Blues",ColAttributes=NULL)
pheatmap(inPeakCors.usingGdna, fontsize=14, main='Using gDNA')
dev.off() # stop putting STDOUT to the file
# It turns out I don't need to worry about pvclust in this case, because it clusters as pheatmap did
#pvlusted.gDNA <- pvclust(inPeakCors.usingGdna, method.hclust='complete')

#plot(pvlusted)
#pdf("pvClustOfCountsInPeaks_subControl.pdf")
#plot(pvlusted.gDNA)
#dev.off()

data.count.contrast.withGDNA = dba.contrast(countsNoL1s, minMembers=2, categories=DBA_CONDITION) # we need to tell DiffBind which fall into which groups
data.count.contrast.analyze.withGDNA = dba.analyze(data.count.contrast.withGDNA, bCorPlot=F, 
                                          bParallel=T, bTagwise=F, bFullLibrarySize=T, 
                                          method=DBA_EDGER, bSubControl=T, bReduceObjects=F) # now finally run the main differential analysis


pdf("PCA_withgDNA_ComBatCorrected.pdf") # this says put the STDOUT to this file
dba.plotPCA(data.count.contrast.analyze.withGDNA, DBA_CONDITION, vColors=c('dodgerblue3' ,'goldenrod3', 'firebrick3')) # the contrast uses the differentially bound sites
dev.off() # stop putting STDOUT to the file

# using that, and the correlation map above, I decided to go with using gDNA

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


pdf("MAPlot_all_smallDots_EDGER_ComBatCorrected.pdf") # this says put the STDOUT to this file
for (i in 1:6){
  dba.plotMA(data.count.contrast.analyze.withGDNA,bSignificant=T,th=qVal,dotSize=0.5, contrast=i, method=DBA_EDGER)
}
dev.off() # stop putting STDOUT to the file


pdf("DifferentialCorrelationHeatmap_all_ComBatCorrected.pdf") # this says put the STDOUT to this file
for (i in 1:6){
  corvals = dba.plotHeatmap(data.count.contrast.analyze.withGDNA, contrast=i, correlations=T)
}
dev.off() # stop putting STDOUT to the file


# Just for reference I want all of the peaks assessed

options(max.print=1e8) 
qVal=1
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
