library(scales)
# All of these number were generated using DiffBind,
#Settings
FC_column=7 # That's 1 based
fdr_column=9
fdrThresh=0.05
setwd('/Users/acd13/Desktop/ATAC/Analysis/differential/N2_dev/allButL1/differentialReports/combatCorrected')

# Input the files
names<-c('Embryo vs. Larval stage 3','Embryo vs. Adult','Larval stage 3 vs. Adult')
files<-c('/Users/acd13/Desktop/ATAC/Analysis/differential/N2_dev/allButL1/differentialReports/combatCorrected/allPeaks/EE_vs_L3_ComBatCorrected_EdgeR_q1.txt',
          '/Users/acd13/Desktop/ATAC/Analysis/differential/N2_dev/allButL1/differentialReports/combatCorrected/allPeaks/EE_vs_YA_ComBatCorrected_EdgeR_q1.txt',
          '/Users/acd13/Desktop/ATAC/Analysis/differential/N2_dev/allButL1/differentialReports/combatCorrected/allPeaks/L3_vs_YA_ComBatCorrected_EdgeR_q1.txt')



file <-read.table(files[1], header=T)
relevantData <-  cbind(-1*file[,FC_column], -1*log10(as.numeric(as.character(file[,fdr_column]))))
pdf('EE_v_L3_volcano_fdr0.05.pdf', width=5, height=5)
plot(relevantData,col=alpha('gray',0.1), pch=20, xlab='Log2(fold change)', ylab='-Log10(FDR)')
points(relevantData[which(relevantData[,2] > -1*log10(fdrThresh) & relevantData[,1] > 0),], col=alpha('goldenrod2',0.3), pch=20)
points(relevantData[which(relevantData[,2] > -1*log10(fdrThresh) & relevantData[,1] < 0),], col=alpha('darkorchid4',0.3), pch=20)
text(5, 40, paste("n=",nrow(relevantData[which(relevantData[,2] >= -1*log10(fdrThresh) & relevantData[,1] > 0),]), sep=""), col='goldenrod2')
text(-5, 40, paste("n=",nrow(relevantData[which(relevantData[,2] >= -1*log10(fdrThresh) & relevantData[,1] < 0),]), sep=""), col='darkorchid4')
text(0, 40, paste("n=",nrow(relevantData[which(relevantData[,2] < -1*log10(fdrThresh)),]), sep=""), col='gray')
dev.off()

file <-read.table(files[2], header=T)
relevantData <-  cbind(-1*file[,FC_column], -1*log10(as.numeric(as.character(file[,fdr_column]))))
pdf('EE_v_YA_volcano_fdr0.05.pdf', width=5, height=5)
plot(relevantData,col=alpha('gray',0.1), pch=20, xlab='Log2(fold change)', ylab='-Log10(FDR)')
points(relevantData[which(relevantData[,2] > -1*log10(fdrThresh) & relevantData[,1] > 0),], col=alpha('darkgreen',0.3), pch=20)
points(relevantData[which(relevantData[,2] > -1*log10(fdrThresh) & relevantData[,1] < 0),], col=alpha('darkorchid4',0.3), pch=20)
text(5, 40, paste("n=",nrow(relevantData[which(relevantData[,2] > -1*log10(fdrThresh) & relevantData[,1] > 0),]), sep=""), col='darkgreen')
text(-5, 40, paste("n=",nrow(relevantData[which(relevantData[,2] > -1*log10(fdrThresh) & relevantData[,1] < 0),]), sep=""), col='darkorchid4')
text(0, 40, paste("n=",nrow(relevantData[which(relevantData[,2] < -1*log10(fdrThresh)),]), sep=""), col='gray')
dev.off()



file <-read.table(files[3], header=T)
relevantData <-  cbind(-1*file[,FC_column], -1*log10(as.numeric(as.character(file[,fdr_column]))))
pdf('L3_v_YA_volcano_fdr0.05.pdf', width=5, height=5)
plot(relevantData,col=alpha('gray',0.1), pch=20, xlab='Log2(fold change)', ylab='-Log10(FDR)')
points(relevantData[which(relevantData[,2] > -1*log10(fdrThresh) & relevantData[,1] > 0),], col=alpha('darkgreen',0.3), pch=20)
points(relevantData[which(relevantData[,2] > -1*log10(fdrThresh) & relevantData[,1] < 0),], col=alpha('goldenrod2',0.3), pch=20)
text(2, 10, paste("n=",nrow(relevantData[which(relevantData[,2] > -1*log10(fdrThresh) & relevantData[,1] > 0),]), sep=""), col='darkgreen')
text(-2, 10, paste("n=",nrow(relevantData[which(relevantData[,2] > -1*log10(fdrThresh) & relevantData[,1] < 0),]), sep=""), col='goldenrod2')
text(0, 10, paste("n=",nrow(relevantData[which(relevantData[,2] < -1*log10(fdrThresh)),]), sep=""), col='gray')
dev.off()

