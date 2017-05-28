library(ggplot2)
library(grid)
#############################
# Define function
#############################
processGorillaAndRevigoFiles <- function(gorillaFile, revigoFile){
  gorillaData <- read.table(gorillaFile, header=T, sep="\t")
  revigoData <- read.csv(revigoFile)
  numberedRevigoData <- numberRevigoGroups(revigoData)
  groupedGorillaData <- merge(gorillaData, numberedRevigoData, by=1, all.x=TRUE)
  
  # format some columns and create some others for plotting
  groupedGorillaData[,'group_number'] <- as.numeric(as.character(groupedGorillaData[,'group_number']))
  groupedGorillaData[,'FDR.q.value'] <- as.numeric(as.character(groupedGorillaData[,'FDR.q.value']))
  groupedGorillaData[,'b'] <- as.numeric(as.character(groupedGorillaData[,'b']))
  groupedGorillaData[,'Enrichment'] <- as.numeric(as.character(groupedGorillaData[,'Enrichment']))
  
  groupedGorillaData$negLog10Fdr <- -log10(groupedGorillaData[,'FDR.q.value'])
  groupedGorillaData$DescriptionWithGeneNum <- paste0(groupedGorillaData$Description, " (", groupedGorillaData$b, ")")
  
  return(groupedGorillaData)
}

numberRevigoGroups <- function(revigoData){
  groupNumber <- 0 # start 1 less than what I actually want the number to be
  toReturn <- cbind(as.character(revigoData[,'term_ID']), rep(0, nrow(revigoData)))
  colnames(toReturn) <- c('term_ID', 'group_number')
  for (i in 1:nrow(revigoData)){
    if(revigoData[i,'eliminated'] == 0){
      groupNumber <- groupNumber + 1
    }
    toReturn[i,'group_number'] <- groupNumber
  }
  return(as.data.frame(toReturn))
}


#############################
# load data
#############################
setwd("/Users/acd13/Desktop/ATAC/Analysis/differential/N2_dev/allButL1/rankingGenesByATACSignalChanges/rankingByAllConnectedAtacPeaks/EE_v_L3/updated2Nov/")
eeGorillaFile <- 'EEvL3_EE_BP_forPlotting.tsv'
eeRevigoFile <- 'EEvL3_EE_BP_forRevigo_results.csv'
eeData <- processGorillaAndRevigoFiles(eeGorillaFile, eeRevigoFile)
eeSorted <- eeData[ order(as.numeric(eeData[,'group_number']), eeData[,'FDR.q.value']), c("DescriptionWithGeneNum", "negLog10Fdr", "Enrichment", "group_number")]
# some of these are still not grouped how I think they should be
# I start from the bottom and add some of the ungrouped or small groups to some of the top groups
# Some that were by themselves and super non-descript I put in 0 
eeSorted[c(89:91), "group_number"] <- 4
eeSorted[which(eeSorted$group_number == 31), "group_number"] <- 4
eeSorted[which(eeSorted$group_number == 30), "group_number"] <- 21
eeSorted[which(eeSorted$group_number == 29), "group_number"] <- 6
eeSorted[which(eeSorted$group_number == 28), "group_number"] <- 1
eeSorted[which(eeSorted$group_number == 27), "group_number"] <- 14
eeSorted[which(eeSorted$group_number == 26), "group_number"] <- 4
eeSorted[which(eeSorted$group_number == 25), "group_number"] <- 1
eeSorted[which(eeSorted$group_number == 24), "group_number"] <- 14
eeSorted[which(eeSorted$group_number == 23), "group_number"] <- 1
eeSorted[which(eeSorted$group_number == 22), "group_number"] <- 4
eeSorted[which(eeSorted$group_number == 20), "group_number"] <- 4
eeSorted[which(eeSorted$group_number == 19), "group_number"] <- 17
eeSorted[which(eeSorted$group_number == 18), "group_number"] <- 6
eeSorted[which(eeSorted$group_number == 15), "group_number"] <- 14
eeSorted[which(eeSorted$group_number == 13), "group_number"] <- 1
eeSorted[which(eeSorted$group_number == 10), "group_number"] <- 11
eeSorted[which(eeSorted$group_number == 9), "group_number"] <- 11
eeSorted[which(eeSorted$group_number == 8), "group_number"] <- 11
eeSorted[which(eeSorted$group_number == 7), "group_number"] <- 1
eeSorted[which(eeSorted$group_number == 5), "group_number"] <- 1
# these were too vague
eeSorted[which(eeSorted$DescriptionWithGeneNum == "biological regulation (28)"), "group_number"] <- NA # original group number 3
eeSorted[which(eeSorted$DescriptionWithGeneNum == "regulation of cellular process (208)"), "group_number"] <- NA # original group number 12
# these only had 1 gene and were in their own group
eeSorted[which(eeSorted$group_number == 16), "group_number"] <- NA

eeResorted <- eeSorted[ order(as.numeric(eeSorted[,'group_number']), -eeSorted[,'negLog10Fdr']), c("DescriptionWithGeneNum", "negLog10Fdr", "Enrichment", "group_number")]
length(unique(eeResorted$group_number))
# 9, which includes 1 that is NA, and we don't want those
# Now I want to group the groups
# The firt super group is development related
devSuperGroupRowInds <- c(which(eeResorted$group_number == 6)
                          ,which(eeResorted$group_number == 21)
                          ,which(eeResorted$group_number == 14)
                          ,which(eeResorted$group_number == 1)
                          )

superGroupOrder <- c(devSuperGroupRowInds
                     ,which(eeResorted$group_number == 4)
                     ,which(eeResorted$group_number == 17)
                     ,which(eeResorted$group_number == 11)
                     ,which(eeResorted$group_number == 2)
                    )

eeFinalSorted <- eeResorted[rev(superGroupOrder),]

# I needed a 2nd pass, but grep wasn't working
rowsToRm <- c(1, 2, 4, 5, 7, 10, 12:14, 16:26, 28:39, 41:44, 46:54, 56,60:65, 67:69, 72:75, 77, 79:83, 85
              , 57 #grep("axon extension (11)", eeFinalSorted$DescriptionWithGeneNum)
              , 27 #grep("transcription, DNA-templated (16)", eeFinalSorted$DescriptionWithGeneNum)
              , 58 #grep("developmental growth involved in morphogenesis (11)", eeFinalSorted$DescriptionWithGeneNum)
              , 55 #grep("developmental cell growth (11)", eeFinalSorted$DescriptionWithGeneNum)
              , 71 #grep("neuron migration (15)", eeFinalSorted$DescriptionWithGeneNum)
              , 6 # behavior
              , 86 # cellular developmental process (108)
) 

write.table(eeFinalSorted[rowsToRm,], file="EE_enriched_termsRemoved.tsv", sep="\t", quote=F, row.names=F)
eeToPlot <- eeFinalSorted[-rowsToRm,]
eeToPlot$DescriptionWithGeneNum <- factor(eeToPlot$DescriptionWithGeneNum, levels = eeToPlot$DescriptionWithGeneNum)

#eeAllPlot = ggplot(eeFinalSorted, aes(x = negLog10Fdr, y = DescriptionWithGeneNum)) + theme_bw(base_size = 18) +
#  xlim(c(0, max(eeFinalSorted$negLog10Fdr)*1.1)) +
#  geom_point(size = log2(eeFinalSorted$Enrichment + 1)*3) +
#  scale_y_discrete(c(1:length(eeFinalSorted$DescriptionWithGeneNum)), labels=eeFinalSorted$DescriptionWithGeneNum) +
#  xlab('-log10(FDR)')

enrichmentScaling <- 4
eeRefinedPlot = ggplot(eeToPlot, aes(x = negLog10Fdr, y = DescriptionWithGeneNum)) + theme_bw(base_size = 18) +
  xlim(c(0, max(eeToPlot$negLog10Fdr)*1.1)) +
  geom_point(size = log2(eeToPlot$Enrichment + 1)*enrichmentScaling) +
  scale_y_discrete(c(1:length(eeToPlot$DescriptionWithGeneNum)), labels=eeToPlot$DescriptionWithGeneNum) +
  xlab('-log10(FDR)')

# Now I plot the scale for
numberDots <- 6
scale.go.term<-c("a1", 'b2',"c3","d4", 'e5', 'f6') # this was a hacky way to get it sorted as I wanted, the letters get removed in Illustator
scale.fold.Enrichment<-rep(1,numberDots)
scale.fdr <- rep(1e-5,numberDots)
scale.geneNumber<-seq(1,numberDots)
scale.negLog10.fdr<- -1*log10(scale.fdr)

scale.df <- data.frame(name = scale.go.term, Enrichment = scale.fold.Enrichment, q = scale.negLog10.fdr, size=scale.geneNumber)
scale.df.sorted.trfrmd <- transform(scale.df, name = factor(name, levels = scale.df$name))

toPlot.df <-data.frame(name = as.character(scale.df.sorted.trfrmd[,1])
                       ,Enrichment=scale.df.sorted.trfrmd$Enrichment
                       , q=scale.df.sorted.trfrmd$q
                       , size=enrichmentScaling*log2(scale.df.sorted.trfrmd$size + 1))                   


# And now combine them
scalePlot = ggplot(toPlot.df, aes(x = q, y = name)) + theme_bw(base_size = 18) +
  xlim(range(toPlot.df$q)) +
  geom_point(size = toPlot.df$size,aes(size= size)) +
  scale_y_discrete(c(1:length(toPlot.df$name)), labels=toPlot.df$name)


# Now we plot everything
# One figure in row 1 and two figures in row 2
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
plot1 <- eeRefinedPlot 
plot2 <- scalePlot

pdf('ggplot2WithScale_EEvL3_rankedByAllConnectedATACPeakTotalChange_SelectGroupedGOProcess_EEenriched.pdf', width=12, height=6)
grid.newpage()
pushViewport(viewport(layout = grid.layout(5, 5)))
print(plot1, vp = vplayout(1:5, 1:4))
print(plot2, vp = vplayout(2:3, 5))
dev.off()






l3GorillaFile <- 'EEvL3_L3_BP_forPlotting.tsv'  
l3RevigoFile <- 'EEvL3_L3_BP_forRevigo_results.csv'
l3Data <- processGorillaAndRevigoFiles(l3GorillaFile, l3RevigoFile)
l3Sorted <- l3Data[ order(as.numeric(l3Data[,'group_number']), l3Data[,'FDR.q.value']), c("DescriptionWithGeneNum", "negLog10Fdr", "Enrichment", "group_number")]
# some of these are still not grouped how I think they should be
# I start from the bottom and add some of the ungrouped or small groups to some of the top groups
# Some that were by themselves and super non-descript I put in 0 
l3Sorted[which(l3Sorted$group_number == 6), "group_number"] <- 5
l3Sorted[which(l3Sorted$group_number == 8), "group_number"] <- 7
l3Sorted[which(l3Sorted$group_number == 10), "group_number"] <- 2
l3Sorted[which(l3Sorted$group_number == 13), "group_number"] <- 1
l3Sorted[which(l3Sorted$group_number == 19), "group_number"] <- 5
l3Sorted[which(l3Sorted$group_number == 20), "group_number"] <- 3
l3Sorted[which(l3Sorted$group_number == 21), "group_number"] <- 24
l3Sorted[which(l3Sorted$group_number == 22), "group_number"] <- 5
l3Sorted[which(l3Sorted$group_number == 23), "group_number"] <- 5
l3Sorted[which(l3Sorted$group_number == 25), "group_number"] <- 5
l3Sorted[which(l3Sorted$group_number == 26), "group_number"] <- 5
l3Sorted[which(l3Sorted$group_number == 27), "group_number"] <- 5
l3Sorted[which(l3Sorted$group_number == 28), "group_number"] <- 5
l3Sorted[which(l3Sorted$group_number == 29), "group_number"] <- 5
l3Sorted[which(l3Sorted$group_number == 31), "group_number"] <- 5
l3Sorted[which(l3Sorted$group_number == 34), "group_number"] <- 5
l3Sorted[which(l3Sorted$group_number == 35), "group_number"] <- 5
l3Sorted[which(l3Sorted$group_number == 36), "group_number"] <- 2
l3Sorted[which(l3Sorted$group_number == 37), "group_number"] <- 5
l3Sorted[which(l3Sorted$group_number == 38), "group_number"] <- 11
l3Sorted[which(l3Sorted$group_number == 40), "group_number"] <- 32
l3Sorted[which(l3Sorted$group_number == 43), "group_number"] <- 5
l3Sorted[which(l3Sorted$group_number == 46), "group_number"] <- 44
l3Sorted[which(l3Sorted$group_number == 47), "group_number"] <- 5
l3Sorted[which(l3Sorted$group_number == 50), "group_number"] <- 5
l3Sorted[which(l3Sorted$group_number == 51), "group_number"] <- 7
l3Sorted[which(l3Sorted$group_number == 52), "group_number"] <- 11
l3Sorted[which(l3Sorted$group_number == 54), "group_number"] <- 44
l3Sorted[which(l3Sorted$group_number == 55), "group_number"] <- 5
l3Sorted[which(l3Sorted$group_number == 56), "group_number"] <- 5
l3Sorted[which(l3Sorted$group_number == 57), "group_number"] <- 5
l3Sorted[which(l3Sorted$group_number == 58), "group_number"] <- 44
l3Sorted[which(l3Sorted$group_number == 59), "group_number"] <- 53
l3Sorted[which(l3Sorted$group_number == 58), "group_number"] <- 7
l3Sorted[which(l3Sorted$group_number == 60), "group_number"] <- 7
l3Sorted[which(l3Sorted$group_number == 61), "group_number"] <- 5
l3Sorted[which(l3Sorted$group_number == 62), "group_number"] <- 5
l3Sorted[which(l3Sorted$group_number == 63), "group_number"] <- 5
l3Sorted[which(l3Sorted$group_number == 64), "group_number"] <- 5
l3Sorted[which(l3Sorted$group_number == 65), "group_number"] <- 5
l3Sorted[which(l3Sorted$group_number == 66), "group_number"] <- 7
l3Sorted[which(l3Sorted$group_number == 67), "group_number"] <- 6
l3Sorted[which(l3Sorted$group_number == 69), "group_number"] <- 11
l3Sorted[which(l3Sorted$group_number == 71), "group_number"] <- 5
l3Sorted[which(l3Sorted$group_number == 72), "group_number"] <- 7
l3Sorted[which(l3Sorted$group_number == 73), "group_number"] <- 39
l3Sorted[which(l3Sorted$group_number == 75), "group_number"] <- 7
l3Sorted[which(l3Sorted$group_number == 76), "group_number"] <- 18
l3Sorted[which(l3Sorted$group_number == 77), "group_number"] <- 7
l3Sorted[which(l3Sorted$group_number == 78), "group_number"] <- 7
l3Sorted[which(l3Sorted$group_number == 79), "group_number"] <- 2
l3Sorted[which(l3Sorted$group_number == 80), "group_number"] <- 2
l3Sorted[which(l3Sorted$group_number == 81), "group_number"] <- 2
l3Sorted[which(l3Sorted$group_number == 82), "group_number"] <- 18
l3Sorted[which(l3Sorted$group_number == 83), "group_number"] <- 39
l3Sorted[which(l3Sorted$group_number == 84), "group_number"] <- 7
l3Sorted[which(l3Sorted$group_number == 85), "group_number"] <- 7
l3Sorted[which(l3Sorted$group_number == 86), "group_number"] <- 2
l3Sorted[which(l3Sorted$group_number == 87), "group_number"] <- 7
l3Sorted[which(l3Sorted$group_number == 88), "group_number"] <- 5
l3Sorted[which(l3Sorted$group_number == 89), "group_number"] <- 1
l3Sorted[which(l3Sorted$group_number == 93), "group_number"] <- 7
l3Sorted[which(l3Sorted$group_number == 95), "group_number"] <- 2
l3Sorted[which(l3Sorted$group_number == 97), "group_number"] <- 5
l3Sorted[which(l3Sorted$group_number == 98), "group_number"] <- 5
l3Sorted[which(l3Sorted$group_number == 99), "group_number"] <- 5
l3Sorted[which(l3Sorted$group_number == 100), "group_number"] <- 1
l3Sorted[which(l3Sorted$group_number == 102), "group_number"] <- 18
l3Sorted[which(l3Sorted$group_number == 103), "group_number"] <- 1
l3Sorted[which(l3Sorted$group_number == 104), "group_number"] <- 14
l3Sorted[which(l3Sorted$group_number == 105), "group_number"] <- 7
l3Sorted[which(l3Sorted$group_number == 106), "group_number"] <- 39
l3Sorted[which(l3Sorted$group_number == 107), "group_number"] <- 7
l3Sorted[which(l3Sorted$group_number == 108), "group_number"] <- 5
l3Sorted[which(l3Sorted$group_number == 109), "group_number"] <- 39
l3Sorted[which(l3Sorted$group_number == 110), "group_number"] <- 18
l3Sorted[which(l3Sorted$group_number == 111), "group_number"] <- 39
l3Sorted[which(l3Sorted$group_number == 112), "group_number"] <- 2
l3Sorted[which(l3Sorted$group_number == 113), "group_number"] <- 2
l3Sorted[which(l3Sorted$group_number == 116), "group_number"] <- 1
l3Sorted[which(l3Sorted$group_number == 117), "group_number"] <- 7
l3Sorted[which(l3Sorted$group_number == 118), "group_number"] <- 3
l3Sorted[which(l3Sorted$group_number == 120), "group_number"] <- 7
l3Sorted[which(l3Sorted$group_number == 122), "group_number"] <- 1
l3Sorted[which(l3Sorted$group_number == 123), "group_number"] <- 7
l3Sorted[which(l3Sorted$group_number == 124), "group_number"] <- 7
l3Sorted[which(l3Sorted$group_number == 127), "group_number"] <- 14
l3Sorted[which(l3Sorted$group_number == 128), "group_number"] <- 2
l3Sorted[which(l3Sorted$group_number == 129), "group_number"] <- 11
l3Sorted[which(l3Sorted$group_number == 130), "group_number"] <- 11
l3Sorted[which(l3Sorted$group_number == 131), "group_number"] <- 3
l3Sorted[which(l3Sorted$group_number == 132), "group_number"] <- 5
l3Sorted[which(l3Sorted$group_number == 133), "group_number"] <- 24
l3Sorted[which(l3Sorted$group_number == 134), "group_number"] <- 2
l3Sorted[which(l3Sorted$group_number == 135), "group_number"] <- 5
l3Sorted[which(l3Sorted$group_number == 138), "group_number"] <- 18
l3Sorted[which(l3Sorted$group_number == 140), "group_number"] <- 1
l3Sorted[which(l3Sorted$group_number == 141), "group_number"] <- 5
l3Sorted[which(l3Sorted$group_number == 142), "group_number"] <- 5
l3Sorted[which(l3Sorted$group_number == 143), "group_number"] <- 44
l3Sorted[which(l3Sorted$group_number == 144), "group_number"] <- 5
l3Sorted[which(l3Sorted$group_number == 145), "group_number"] <- 5
l3Sorted[which(l3Sorted$group_number == 146), "group_number"] <- 2
l3Sorted[which(l3Sorted$group_number == 147), "group_number"] <- 2
l3Sorted[which(l3Sorted$group_number == 148), "group_number"] <- 5
l3Sorted[which(l3Sorted$group_number == 149), "group_number"] <- 44
l3Sorted[which(l3Sorted$group_number == 150), "group_number"] <- 5
l3Sorted[which(l3Sorted$group_number == 152), "group_number"] <- 1
l3Sorted[which(l3Sorted$group_number == 156), "group_number"] <- 18
l3Sorted[which(l3Sorted$group_number == 157), "group_number"] <- 18
l3Sorted[which(l3Sorted$group_number == 158), "group_number"] <- 18
l3Sorted[which(l3Sorted$group_number == 159), "group_number"] <- 1
l3Sorted[which(l3Sorted$group_number == 161), "group_number"] <- 18
l3Sorted[which(l3Sorted$group_number == 160), "group_number"] <- 7
l3Sorted[which(l3Sorted$group_number == 162), "group_number"] <- 5

# these were too vague
l3Sorted[which(l3Sorted$DescriptionWithGeneNum == "cellular process (390)"), "group_number"] <- NA # original group number 6
l3Sorted[which(l3Sorted$DescriptionWithGeneNum == "multicellular organismal process (105)"), "group_number"] <- NA # original group number 9
l3Sorted[which(l3Sorted$DescriptionWithGeneNum == "single-organism process (2409)"), "group_number"] <- NA # original group number 12
# I got tired of copying and pasting the vague names
l3Sorted[which(l3Sorted$group_number == 15), "group_number"] <- NA
l3Sorted[which(l3Sorted$group_number == 16), "group_number"] <- NA
l3Sorted[which(l3Sorted$group_number == 17), "group_number"] <- NA
l3Sorted[which(l3Sorted$group_number == 33), "group_number"] <- NA
l3Sorted[which(l3Sorted$group_number == 115), "group_number"] <- NA
l3Sorted[which(l3Sorted$group_number == 4), "group_number"] <- NA # behavior
l3Sorted[which(l3Sorted$group_number == 14), "group_number"] <- NA # localization
l3Sorted[which(l3Sorted$group_number == 45), "group_number"] <- NA #negative regulation of biological process
l3Sorted[which(l3Sorted$group_number == 48), "group_number"] <- NA #positive regulation of biological process
l3Sorted[which(l3Sorted$group_number == 68), "group_number"] <- NA #regulation of cellular process
l3Sorted[which(l3Sorted$group_number == 91), "group_number"] <- NA # system process
l3Sorted[which(l3Sorted$group_number == 121), "group_number"] <- NA # single-multicellular organism process
l3Sorted[which(l3Sorted$group_number == 125), "group_number"] <- NA # regulation of biological process


l3Resorted <- l3Sorted[ order(as.numeric(l3Sorted[,'group_number']), -l3Sorted[,'negLog10Fdr']), c("DescriptionWithGeneNum", "negLog10Fdr", "Enrichment", "group_number")]
length(unique(l3Resorted$group_number))
# 35 which includes 1 that is NA, and we don't want those
# Now I want to group the groups
# The firt super group is larval development related
larvDevSuperGroupRowInds <- c(which(l3Resorted$group_number == 2)
                          ,which(l3Resorted$group_number == 1)
                          ,which(l3Resorted$group_number == 11)
                          ,which(l3Resorted$group_number == 24)
)

# The second super group is signaling related
sigSuperGroupRowInds <- c(which(l3Resorted$group_number == 7)
                              ,which(l3Resorted$group_number == 44)
                              ,which(l3Resorted$group_number == 53)
)

# The third super group is metabolism related
metabSuperGroupRowInds <- c(which(l3Resorted$group_number == 5)
                          ,which(l3Resorted$group_number == 18)
)

# The last super group is quality control and aging
#sigSuperGroupRowInds <- c(which(l3Resorted$group_number == 32) # qc
#                          ,which(l3Resorted$group_number == 101) # aging
#)
write.table(l3Resorted, file="L3_enriched_allGroupedTerms.tsv", sep="\t", quote=F, row.names=F)

superGroupOrder <- c(larvDevSuperGroupRowInds
                     ,sigSuperGroupRowInds
                     ,metabSuperGroupRowInds
                     ,which(l3Resorted$group_number == 39) # gene expression reg
                     ,which(l3Resorted$group_number == 32) # qc
                     ,which(l3Resorted$group_number == 101) # aging
)

l3FinalSorted <- l3Resorted[rev(superGroupOrder),]

# I also decided to switch up the order after all this...
#rowsToPlot <- c(1, 11, 12, 14, 42, 43, 84, 92, 119, 126, 131, 139, 151, 153, 162, 175, 179, 183, 201, 208, 215, 221, 229) 
rowsToPlot <- c(1, 12, 42, 84, 92, 119, 126, 151, 131, 139, 162, 179, 183, 201, 208, 215, 175, 221, 229) 

l3ToPlot <- l3FinalSorted[rowsToPlot,]
l3ToPlot$DescriptionWithGeneNum <- factor(l3ToPlot$DescriptionWithGeneNum, levels = l3ToPlot$DescriptionWithGeneNum)

#l3AllPlot = ggplot(l3FinalSorted, aes(x = negLog10Fdr, y = DescriptionWithGeneNum)) + theme_bw(base_size = 18) +
#  xlim(c(0, max(l3FinalSorted$negLog10Fdr)*1.1)) +
#  geom_point(size = log2(l3FinalSorted$Enrichment + 1)*3) +
#  scale_y_discrete(c(1:length(l3FinalSorted$DescriptionWithGeneNum)), labels=l3FinalSorted$DescriptionWithGeneNum) +
#  xlab('-log10(FDR)')

enrichmentScaling <- 4
l3RefinedPlot = ggplot(l3ToPlot, aes(x = negLog10Fdr, y = DescriptionWithGeneNum)) + theme_bw(base_size = 18) +
  xlim(c(0, max(l3ToPlot$negLog10Fdr)*1.1)) +
  geom_point(size = log2(l3ToPlot$Enrichment + 1)*enrichmentScaling) +
  scale_y_discrete(c(1:length(l3ToPlot$DescriptionWithGeneNum)), labels=l3ToPlot$DescriptionWithGeneNum) +
  xlab('-log10(FDR)')

# Now I plot the scale for
numberDots <- 6
scale.go.term<-c("a1", 'b2',"c3","d4", 'e5', 'f6') # this was a hacky way to get it sorted as I wanted, the letters get removed in Illustator
scale.fold.Enrichment<-rep(1,numberDots)
scale.fdr <- rep(1e-5,numberDots)
scale.geneNumber<-seq(1,numberDots)
scale.negLog10.fdr<- -1*log10(scale.fdr)

scale.df <- data.frame(name = scale.go.term, Enrichment = scale.fold.Enrichment, q = scale.negLog10.fdr, size=scale.geneNumber)
scale.df.sorted.trfrmd <- transform(scale.df, name = factor(name, levels = scale.df$name))

toPlot.df <-data.frame(name = as.character(scale.df.sorted.trfrmd[,1])
                       ,Enrichment=scale.df.sorted.trfrmd$Enrichment
                       , q=scale.df.sorted.trfrmd$q
                       , size=enrichmentScaling*log2(scale.df.sorted.trfrmd$size + 1))                   


# And now combine them
scalePlot = ggplot(toPlot.df, aes(x = q, y = name)) + theme_bw(base_size = 18) +
  xlim(range(toPlot.df$q)) +
  geom_point(size = toPlot.df$size,aes(size= size)) +
  scale_y_discrete(c(1:length(toPlot.df$name)), labels=toPlot.df$name)


# Now we plot everything
# One figure in row 1 and two figures in row 2
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
plot1 <- l3RefinedPlot 
plot2 <- scalePlot

pdf('ggplot2WithScale_EEvL3_rankedByAllConnectedATACPeakTotalChange_SelectGroupedGOProcess_L3enriched.pdf', width=12, height=6)
grid.newpage()
pushViewport(viewport(layout = grid.layout(5, 5)))
print(plot1, vp = vplayout(1:5, 1:4))
print(plot2, vp = vplayout(2:3, 5))
dev.off()
