library(fields)
library(dplyr)

#############
# Function
fudgeit <- function(){
  xm <- get('xm', envir = parent.frame(1))
  ym <- get('ym', envir = parent.frame(1))
  z  <- get('dens', envir = parent.frame(1))
  colramp <- get('colramp', parent.frame(1))
  image.plot(xm,ym,z, col = colramp(256), legend.only = T, add =F)
}

Break2ColsIntoPercentiles <- function (matrix, colOfInterest){
  binNum<-20
  bins<-seq(0,1,(1/binNum))
  
  colOfInterestPercentiles<-vector(mode='list', length=binNum)
  otherColPercentiles<-vector(mode='list', length=binNum)
  
  for (i in 1:binNum){
    
    low <- quantile(matrix[,colOfInterest], bins[i])
    high <-quantile(matrix[,colOfInterest], bins[i+1])
    
    colOfInterestPercentiles[[i]]<- matrix[which(matrix[,colOfInterest]<=high & matrix[,colOfInterest]>=low),colOfInterest]
    otherColPercentiles[[i]]<- matrix[which(matrix[,colOfInterest]<=high & matrix[,colOfInterest]>=low),(-1*colOfInterest)]
    
  }
  toReturn<-vector(mode='list', length=2)
  toReturn$colOfInterest=colOfInterestPercentiles
  toReturn$otherCol=otherColPercentiles
  
  return(toReturn)
  
}


####################
setwd("/Users/acd13/Desktop/ATAC/Analysis/gro/expression/")

# read in the data, being sure to sort by gene name (which when I generated the coverage data I amde sure was unique and descriptive)
# First the gro-seq coverage
l3_GROCoverage_refseq <- read.csv('/Users/acd13/Desktop/GROseq/L3/quantification/L3Merged.homerAnalyzeRepeats_refseqRNA_strandSpecific_ce10_downstreamReadThrough1kb.processed.txt', sep='\t', header=T)
# And eliminate duplicate gene names (presumably from isoforms)
l3_GROCoverage_refseq <- l3_GROCoverage_refseq %>%  group_by(Gene) %>%  summarise(max(FPKM, na.rm=TRUE))

# Then the atac-seq promoter 
allATACCoverage_refChenExtendedProms <- read.csv('/Users/acd13/Desktop/ATAC/Analysis/coverage/inGenicLoci/inProms/EE_L1_L3_YA_gDNA_merged_inserts_InExtProms.bed', sep='\t', header=F)
colnames(allATACCoverage_refChenExtendedProms) <- c("chr","start","stop","geneName","ignore","strand","EE","L1","L3","YA","gDNA")

# And eliminate duplicate gene names (presumably from isoforms)
l3ATACCoverage_refChenExtendedProms <- allATACCoverage_refChenExtendedProms %>%  group_by(geneName) %>%  summarise(max(L3, na.rm=TRUE))

# combine things for ease
L3GroAndExtAtac <- merge(l3_GROCoverage_refseq, l3ATACCoverage_refChenExtendedProms, by=1)

# Split into percentiles by ATAC-seq signal in promoter
L3GroReadThroughAndExtATAC_PercentilesByPromsATAC<-Break2ColsIntoPercentiles(L3GroAndExtAtac[,c(2,3)],2)

percentileNames<-c("0-5th",'5-10th','10-15th','15-20th','20-25th','25-30th','30-35th','35-40th','40-45th','45-50th','50-55th','55-60th','60-65th','65-70th','70-75th','75-80th','80-85th','85-90th','90-95th','95-100th')

pdf('ExtPromL3ATACSignal_and_L3GRO_signal_refChenGenes_splitIntoPercentilesByATACSignalInProm.pdf')
  boxplot(L3GroReadThroughAndExtATAC_PercentilesByPromsATAC$colOfInterest, notch=T, outline=F,names=percentileNames, ylab='Total ATAC inserts in 1.5kb promoter', xlab='Promoter ATAC-signal percentiles', col='goldenrod3')
  boxplot(L3GroReadThroughAndExtATAC_PercentilesByPromsATAC$otherCol, notch=T, outline=F,names=percentileNames, ylab='GRO-seq FPKM', xlab='1.5kb Promoter ATAC-signal percentiles', col='goldenrod3')
dev.off()

# I also want the scatter plots, just in case
# L3 next
toPlot<-cbind(log10(L3GroAndExtAtac[,2]+0.01), L3GroAndExtAtac[,3])

correlation<-round(cor(L3GroAndExtAtac[,2], L3GroAndExtAtac[,3]), digits=3)
textToAdd<-paste('R= ',correlation,sep='')

xmax<-max(toPlot[,1])
xmin<-max(min(toPlot[,1]),0)

ymax<-max(toPlot[,2])
ymin<-max(min(toPlot[,2]),0)

xlab.name="Log10(total ATAC inserts/1.5kb promoter)"
ylab.name="GRO-seq FPKM"
Lab.palette <- colorRampPalette(c("white", "goldenrod3", 'black'), space = "Lab")

pdf("L3GROseq_as_function_of_L3ATAC_signal_in1.5kbProm_refSeqGenes.pdf")
smoothScatter(toPlot, colramp = Lab.palette, xlab = xlab.name, ylab = ylab.name, 
              main = "", xlim = c(xmin,xmax),ylim = c(ymin,ymax),col="grey", cex=1, pch=20,postPlotHook = fudgeit)
legend('topleft', textToAdd, bty='n')
dev.off()
