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
setwd("/Users/acd13/Desktop/ATAC/Analysis/differential/N2_dev/allButL1/rankingGenesByATACSignalChanges/rankingByAllConnectedAtacPeaks/L3_v_YA/updated2Nov/")
l3GorillaFile <- 'L3vYA_L3_BP_forPlotting.tsv'
l3RevigoFile <- 'L3vYA_L3_BP_forRevigo_results.csv'
l3Data <- processGorillaAndRevigoFiles(l3GorillaFile, l3RevigoFile)
l3Sorted <- l3Data[ order(as.numeric(l3Data[,'group_number']), l3Data[,'FDR.q.value']), c("DescriptionWithGeneNum", "negLog10Fdr", "Enrichment", "group_number")]
# some of these are still not grouped how I think they should be
# I start from the bottom and add some of the ungrouped or small groups to some of the top groups
# Some that were by themselves and super non-descript I put in 0 
l3Sorted[which(is.na(l3Sorted$group_number)), "group_number"] <- 84 # there were three all related to transcription
l3Sorted[which(l3Sorted$group_number == 29), "group_number"] <- 84
l3Sorted[which(l3Sorted$group_number == 71), "group_number"] <- 84
l3Sorted[which(l3Sorted$group_number == 66), "group_number"] <- 84

l3Sorted[which(l3Sorted$group_number == 3), "group_number"] <- 4 #reproduction
l3Sorted[which(l3Sorted$group_number == 32), "group_number"] <- 4 

l3Sorted[which(l3Sorted$group_number == 39), "group_number"] <- 38 # molting

l3Sorted[which(l3Sorted$group_number == 24), "group_number"] <- 62 # cell homeostasis/qc

l3Sorted[which(l3Sorted$group_number == 109), "group_number"] <- 33 # stimulus response
l3Sorted[which(l3Sorted$group_number == 36), "group_number"] <- 33 #
l3Sorted[which(l3Sorted$group_number == 107), "group_number"] <- 33 # 

l3Sorted[which(l3Sorted$group_number == 100), "group_number"] <- 34 # locomotion/muscle development
l3Sorted[which(l3Sorted$group_number == 102), "group_number"] <- 34
l3Sorted[which(l3Sorted$group_number == 108), "group_number"] <- 34
l3Sorted[which(l3Sorted$group_number == 96), "group_number"] <- 34
l3Sorted[which(l3Sorted$group_number == 34), "group_number"] <- 9
l3Sorted[which(l3Sorted$group_number == 30), "group_number"] <- 9
l3Sorted[which(l3Sorted$group_number == 51), "group_number"] <- 9
l3Sorted[which(l3Sorted$group_number == 68), "group_number"] <- 9
l3Sorted[which(l3Sorted$group_number == 81), "group_number"] <- 9
l3Sorted[which(l3Sorted$group_number == 87), "group_number"] <- 9

l3Sorted[which(l3Sorted$group_number == 106), "group_number"] <- 1 # cell component movement
l3Sorted[which(l3Sorted$group_number == 25), "group_number"] <- 1
l3Sorted[which(l3Sorted$group_number == 86), "group_number"] <- 1
l3Sorted[which(l3Sorted$group_number == 91), "group_number"] <- 1
l3Sorted[which(l3Sorted$group_number == 101), "group_number"] <- 1
l3Sorted[which(l3Sorted$group_number == 61), "group_number"] <- 1
l3Sorted[which(l3Sorted$group_number == 14), "group_number"] <- 1

l3Sorted[which(l3Sorted$group_number == 5), "group_number"] <- 2 # cell adhesion
l3Sorted[which(l3Sorted$group_number == 15), "group_number"] <- 2

l3Sorted[which(l3Sorted$group_number == 18), "group_number"] <- 6 # cell communication and signaling
l3Sorted[which(l3Sorted$group_number == 31), "group_number"] <- 6
l3Sorted[which(l3Sorted$group_number == 35), "group_number"] <- 6
l3Sorted[which(l3Sorted$group_number == 40), "group_number"] <- 6
l3Sorted[which(l3Sorted$group_number == 46), "group_number"] <- 6
l3Sorted[which(l3Sorted$group_number == 49), "group_number"] <- 6
l3Sorted[which(l3Sorted$group_number == 50), "group_number"] <- 6
l3Sorted[which(l3Sorted$group_number == 54), "group_number"] <- 6
l3Sorted[which(l3Sorted$group_number == 57), "group_number"] <- 6
l3Sorted[which(l3Sorted$group_number == 58), "group_number"] <- 6
l3Sorted[which(l3Sorted$group_number == 64), "group_number"] <- 6
l3Sorted[which(l3Sorted$group_number == 67), "group_number"] <- 6
l3Sorted[which(l3Sorted$group_number == 69), "group_number"] <- 6
l3Sorted[which(l3Sorted$group_number == 75), "group_number"] <- 6
l3Sorted[which(l3Sorted$group_number == 76), "group_number"] <- 6
l3Sorted[which(l3Sorted$group_number == 82), "group_number"] <- 6
l3Sorted[which(l3Sorted$group_number == 44), "group_number"] <- 6
l3Sorted[which(l3Sorted$group_number == 80), "group_number"] <- 6

l3Sorted[which(l3Sorted$group_number == 37), "group_number"] <- 72 # dauer

l3Sorted[which(l3Sorted$group_number == 59), "group_number"] <- 94 # biosynthesis

l3Sorted[which(l3Sorted$group_number == 60), "group_number"] <- 105 # metabloism regulation

l3Sorted[which(l3Sorted$group_number == 95), "group_number"] <- 45# growth and developemtn and morphogenesis
l3Sorted[which(l3Sorted$group_number == 63), "group_number"] <- 45
l3Sorted[which(l3Sorted$group_number == 65), "group_number"] <- 45
l3Sorted[which(l3Sorted$group_number == 56), "group_number"] <- 45
l3Sorted[which(l3Sorted$group_number == 90), "group_number"] <- 45
l3Sorted[which(l3Sorted$group_number == 23), "group_number"] <- 45 
l3Sorted[which(l3Sorted$group_number == 22), "group_number"] <- 45 
l3Sorted[which(l3Sorted$group_number == 8), "group_number"] <- 45
l3Sorted[which(l3Sorted$group_number == 73), "group_number"] <- 45
l3Sorted[which(l3Sorted$group_number == 74), "group_number"] <- 45
l3Sorted[which(l3Sorted$group_number == 78), "group_number"] <- 45
l3Sorted[which(l3Sorted$group_number == 97), "group_number"] <- 45
l3Sorted[which(l3Sorted$group_number == 99), "group_number"] <- 45
l3Sorted[which(l3Sorted$group_number == 104), "group_number"] <- 45
l3Sorted[which(l3Sorted$group_number == 93), "group_number"] <- 45
l3Sorted[which(l3Sorted$group_number == 70), "group_number"] <- 45
l3Sorted[which(l3Sorted$group_number == 20), "group_number"] <- 45
l3Sorted[which(l3Sorted$group_number == 77), "group_number"] <- 45

l3Sorted[which(l3Sorted$group_number == 17), "group_number"] <- 89 # actin related
l3Sorted[which(l3Sorted$group_number == 55), "group_number"] <- 89
l3Sorted[which(l3Sorted$group_number == 88), "group_number"] <- 89

l3Sorted[which(l3Sorted$group_number == 52), "group_number"] <- 47 # axon/cell extension/projection
l3Sorted[which(l3Sorted$group_number == 103), "group_number"] <- 47
l3Sorted[which(l3Sorted$group_number == 92), "group_number"] <- 47
l3Sorted[which(l3Sorted$group_number == 98), "group_number"] <- 47
l3Sorted[which(l3Sorted$group_number == 85), "group_number"] <- 47
l3Sorted[which(l3Sorted$group_number == 43), "group_number"] <- 47

# these were too vague
l3Sorted[which(l3Sorted$DescriptionWithGeneNum == "multicellular organismal process (94)"), "group_number"] <- NA 
l3Sorted[which(l3Sorted$DescriptionWithGeneNum == "single-organism process (109)"), "group_number"] <- NA 
l3Sorted[which(l3Sorted$DescriptionWithGeneNum == "regulation of multicellular organismal process (123)"), "group_number"] <- NA 
l3Sorted[which(l3Sorted$DescriptionWithGeneNum == "multi-organism process (90)"), "group_number"] <- NA 
l3Sorted[which(l3Sorted$DescriptionWithGeneNum == "biological regulation (475)"), "group_number"] <- NA 
l3Sorted[which(l3Sorted$DescriptionWithGeneNum == "single-organism cellular process (50)"), "group_number"] <- NA 
l3Sorted[which(l3Sorted$DescriptionWithGeneNum == "negative regulation of biological process (130)"), "group_number"] <- NA 
l3Sorted[which(l3Sorted$DescriptionWithGeneNum == "regulation of molecular function (122)"), "group_number"] <- NA 
l3Sorted[which(l3Sorted$DescriptionWithGeneNum == "regulation of cellular process (361)"), "group_number"] <- NA 
l3Sorted[which(l3Sorted$DescriptionWithGeneNum == "positive regulation of biological process (189)"), "group_number"] <- NA 
l3Sorted[which(l3Sorted$DescriptionWithGeneNum == "system process (54)"), "group_number"] <- NA 
l3Sorted[which(l3Sorted$DescriptionWithGeneNum == "regulation of biological process (261)"), "group_number"] <- NA 

l3Resorted <- l3Sorted[ order(as.numeric(l3Sorted[,'group_number']), -l3Sorted[,'negLog10Fdr']), c("DescriptionWithGeneNum", "negLog10Fdr", "Enrichment", "group_number")]
length(unique(l3Resorted$group_number))
# 21, which includes 1 that is NA, and we don't want those
# Now I want to group the groups
# The firt super group is development related
devSuperGroupRowInds <- c(which(l3Resorted$group_number == 72) # dauer
                          ,which(l3Resorted$group_number == 45) # morphogenesis/dev
                          ,which(l3Resorted$group_number == 47) # nerv sys dev
                          ,which(l3Resorted$group_number == 38) # molting
                          ,which(l3Resorted$group_number == 21) # cell death
)
cellProcess <- c(which(l3Resorted$group_number == 1) # intracel movement
                          ,which(l3Resorted$group_number == 2) # cell adhesion
                          ,which(l3Resorted$group_number == 89) # actin related
                          ,which(l3Resorted$group_number == 16) # polarity
)
superGroupOrder <- c(devSuperGroupRowInds
                     ,which(l3Resorted$group_number == 6)
                     ,which(l3Resorted$group_number == 4)
                     ,which(l3Resorted$group_number == 9)
                     ,cellProcess
                    )

l3FinalSorted <- l3Resorted[rev(superGroupOrder),]

# I just pulled out those of interest, and re-ordered for consistency of signal decrease
toPlot <- c(60,22, 28, 37,47,48, 51, 74,83,94,97,143,126,173) 

write.table(l3FinalSorted[-toPlot,], file="l3_enriched_termsRemoved.tsv", sep="\t", quote=F, row.names=F)
l3ToPlot <- l3FinalSorted[toPlot,]
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
numberDots <- 3
scale.go.term<-c('b2',"d4", 'f6') # this was a hacky way to get it sorted as I wanted, the letters get removed in Illustator
scale.fold.Enrichment<-rep(1,numberDots)
scale.fdr <- rep(1e-5,numberDots)
scale.geneNumber<-seq(2,numberDots*2, 2)
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

pdf('ggplot2WithScale_L3vYA_rankedByAllConnectedATACPeakTotalChange_SelectGroupedGOProcess_L3enriched.pdf', width=12, height=6)
grid.newpage()
pushViewport(viewport(layout = grid.layout(5, 5)))
print(plot1, vp = vplayout(1:5, 1:4))
print(plot2, vp = vplayout(2:3, 5))
dev.off()






yaGorillaFile <- 'L3vYA_YA_BP_forPlotting.tsv'  
yaRevigoFile <- 'L3vYA_YA_BP_forRevigo_results.csv'
yaData <- processGorillaAndRevigoFiles(yaGorillaFile, yaRevigoFile)
yaSorted <- yaData[ order(as.numeric(yaData[,'group_number']), yaData[,'FDR.q.value']), c("DescriptionWithGeneNum", "negLog10Fdr", "Enrichment", "group_number")]

# some of these are still not grouped how I think they should be
# I start from the bottom and add some of the ungrouped or small groups to some of the top groups
# Some that were by themselves and super non-descript I put in 0 

yaSorted[which(is.na(yaSorted$group_number)), "group_number"] <- 11 # the only one was DNA segragation
yaSorted[which(yaSorted$group_number == 12), "group_number"] <- 11 # meiosis/cell cycle
yaSorted[which(yaSorted$group_number == 14), "group_number"] <- 11
yaSorted[which(yaSorted$group_number == 8), "group_number"] <- 11
yaSorted[which(yaSorted$group_number == 28), "group_number"] <- 11
yaSorted[which(yaSorted$group_number == 20), "group_number"] <- 11
yaSorted[which(yaSorted$group_number == 33), "group_number"] <- 11
yaSorted[which(yaSorted$group_number == 35), "group_number"] <- 11
yaSorted[which(yaSorted$group_number == 43), "group_number"] <- 11
yaSorted[which(yaSorted$group_number == 41), "group_number"] <- 11

yaSorted[which(yaSorted$group_number == 4), "group_number"] <- 1 # reprod
yaSorted[which(yaSorted$group_number == 32), "group_number"] <- 1 
yaSorted[which(yaSorted$group_number == 34), "group_number"] <- 1 
yaSorted[which(yaSorted$group_number == 36), "group_number"] <- 1 
yaSorted[which(yaSorted$group_number == 9), "group_number"] <- 1 

yaSorted[which(yaSorted$group_number == 17), "group_number"] <- 16 # metab
yaSorted[which(yaSorted$group_number == 18), "group_number"] <- 16 # 
yaSorted[which(yaSorted$group_number == 15), "group_number"] <- 16 # 

yaSorted[which(yaSorted$group_number == 2), "group_number"] <- 26 # rna processing/metabloism
yaSorted[which(yaSorted$group_number == 30), "group_number"] <- 26 
yaSorted[which(yaSorted$group_number == 31), "group_number"] <- 26 
yaSorted[which(yaSorted$group_number == 24), "group_number"] <- 26 

yaSorted[which(yaSorted$group_number == 27), "group_number"] <- 23 # postranscriptional silencing

yaSorted[which(yaSorted$group_number == 29), "group_number"] <- 21 # translation

yaSorted[which(yaSorted$group_number == 38), "group_number"] <- 37 # dev
yaSorted[which(yaSorted$group_number == 40), "group_number"] <- 37
yaSorted[which(yaSorted$group_number == 3), "group_number"] <- 37
yaSorted[which(yaSorted$group_number == 5), "group_number"] <- 37
yaSorted[which(yaSorted$group_number == 39), "group_number"] <- 37
yaSorted[which(yaSorted$group_number == 42), "group_number"] <- 37
yaSorted[which(yaSorted$group_number == 25), "group_number"] <- 37
yaSorted[which(yaSorted$group_number == 10), "group_number"] <- 37
yaSorted[which(yaSorted$group_number == 13), "group_number"] <- 37

yaSorted[which(yaSorted$group_number == 19), "group_number"] <- 7 # dna/general stress response
yaSorted[which(yaSorted$group_number == 6), "group_number"] <- 7 

yaResorted <- yaSorted[ order(as.numeric(yaSorted[,'group_number']), -yaSorted[,'negLog10Fdr']), c("DescriptionWithGeneNum", "negLog10Fdr", "Enrichment", "group_number")]
length(unique(yaResorted$group_number))
# 11 which includes 1 that is NA, and we don't want those

write.table(yaResorted, file="ya_enriched_allGroupedTerms.tsv", sep="\t", quote=F, row.names=F)

superGroupOrder <- c(which(yaResorted$group_number == 1) # Reprod
                     ,which(yaResorted$group_number == 11) # meiosis
                     ,which(yaResorted$group_number == 7) # dna stress response
                     ,which(yaResorted$group_number == 37) # dev
                     ,which(yaResorted$group_number == 26) # rna processing/metabloism
                     ,which(yaResorted$group_number == 22) # gene regulation
                     ,which(yaResorted$group_number == 23) # postranscriptional silencing
                     ,which(yaResorted$group_number == 21) # translation
                     ,which(yaResorted$group_number == 16) # metab
)

yaFinalSorted <- yaResorted[rev(superGroupOrder),]

# I also decided to switch up the order after all this...
rowsToPlot <- c(5, 8,11, 20, 21, 27, 33, 34, 53, 57, 59, 63, 67, 70) 

yaToPlot <- yaFinalSorted[rowsToPlot,]
yaToPlot$DescriptionWithGeneNum <- factor(yaToPlot$DescriptionWithGeneNum, levels = yaToPlot$DescriptionWithGeneNum)

#yaAllPlot = ggplot(yaFinalSorted, aes(x = negLog10Fdr, y = DescriptionWithGeneNum)) + theme_bw(base_size = 18) +
#  xlim(c(0, max(yaFinalSorted$negLog10Fdr)*1.1)) +
#  geom_point(size = log2(yaFinalSorted$Enrichment + 1)*3) +
#  scale_y_discrete(c(1:length(yaFinalSorted$DescriptionWithGeneNum)), labels=yaFinalSorted$DescriptionWithGeneNum) +
#  xlab('-log10(FDR)')

enrichmentScaling <- 4
yaRefinedPlot = ggplot(yaToPlot, aes(x = negLog10Fdr, y = DescriptionWithGeneNum)) + theme_bw(base_size = 18) +
  xlim(c(0, max(yaToPlot$negLog10Fdr)*1.1)) +
  geom_point(size = log2(yaToPlot$Enrichment + 1)*enrichmentScaling) +
  scale_y_discrete(c(1:length(yaToPlot$DescriptionWithGeneNum)), labels=yaToPlot$DescriptionWithGeneNum) +
  xlab('-log10(FDR)')

# Now I plot the scale for
numberDots <- 3
scale.go.term<-c('b2',"d4", 'f6') # this was a hacky way to get it sorted as I wanted, the letters get removed in Illustator
scale.fold.Enrichment<-rep(1,numberDots)
scale.fdr <- rep(1e-5,numberDots)
scale.geneNumber<-seq(2,numberDots*2, 2)
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
plot1 <- yaRefinedPlot 
plot2 <- scalePlot

pdf('ggplot2WithScale_L3vYA_rankedByAllConnectedATACPeakTotalChange_SelectGroupedGOProcess_YAenriched.pdf', width=12, height=6)
grid.newpage()
pushViewport(viewport(layout = grid.layout(5, 5)))
print(plot1, vp = vplayout(1:5, 1:4))
print(plot2, vp = vplayout(2:3, 5))
dev.off()
