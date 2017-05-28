###### Functions ##########
calc.mean.sd_as_error <- function(my.table) {
  rowavg <- rowMeans(my.table)
  sd <-apply(my.table,1,sd)
  upper <- (rowavg+sd)
  lower <- (rowavg-sd)
  my.ouput <- cbind(rowavg,upper,lower)
}

parseResults <- function(fileNames,comparNames, SampleNames){
  nComps<- length(comparNames)
  
  # matrices for later
  Percs <- pvals <- MedianPercs <- maxPercs <- minPercs <- Enrich <- matrix(nrow=nComps,ncol=length(SampleNames),byrow=T) #create an empty matrix to dump results into
  rownames(Percs) <- rownames(pvals) <- rownames(minPercs) <- rownames(MedianPercs) <- rownames(maxPercs) <- rownames(Enrich)<- comparNames
  colnames(Percs) <- colnames(pvals) <- colnames(minPercs) <- colnames(MedianPercs) <- colnames(maxPercs) <- colnames(Enrich)<- SampleNames
  
  for(i in 1:length(SampleNames)){
    for(j in 1:length(comparNames)){ # needed because the files contain replicates and I want to go through them at the same time
      
      acutalFileNumber<-i+((j-1)*length(SampleNames)) # because of the replicates, and they're in order: all of comparisons with rep1, all comparisons with rep2...
      
      featureName<-comparNames[i]
      my.file <- fileNames[acutalFileNumber]
      
      # set the title
      sampleName<-SampleNames[j]
      title<-paste(sampleName,featureName,sep="_")
      
      # read in file which contains the null distribution, make sure it is numeric
      # I wrote a perl script which uses bedtools to do all of this
      NullDist <- scan(my.file, what = 'numeric', sep="\n")
      total <- as.numeric(substring(NullDist[1],2)) # the first entry is always the number of total peaks, but has a '#' in front of it
      actual <-  as.numeric(substring(NullDist[2],2)) # the 2nd entry is the actual number of peaks, but has a '#' in front of it
      
      # get the null dist in the form I want
      NullDist<-NullDist[-2] # this removes the first entry, which is a header containing 'actual'
      NullDist <- as.numeric(NullDist)
      # now do the math
      actualPerc <- (actual/total)*100
      actualPerc <- signif(actualPerc, 4)
      Percs[j,i]<-actualPerc
      
      NullDistPerc <- (NullDist/total)*100
      NullDistPercMax <- signif(max(NullDistPerc, na.rm=T), 4)
      NullDistPercMin <- signif(min(NullDistPerc, na.rm=T), 4)
      NullDistPercMedian <- signif(median(NullDistPerc, na.rm=T), 4)
      MedianPercs[j,i]<-NullDistPercMedian
      maxPercs[j,i]<-NullDistPercMax
      minPercs[j,i]<-NullDistPercMin
      
      # calculate the Fold Enrichment
      if(NullDistPercMedian==0) fc <- 0.1 else fc <- signif(actualPerc/NullDistPercMedian,3)
      
      Enrich[j,i]<-fc
      
      # Calculate the p-value, if 0, it means it is less than 1/(number of times shuffled (from perl file))
      Fn <- ecdf(NullDist)
      if (log2(fc)>0){ # if it's enriched we're on the right tail, so it's 1-
        pval <- 1- Fn(actual)
      }else{ # otherwise it's depleted, so we just take the valu
        pval <- Fn(actual)
      }
      
      pvals[j,i]<-pval
    }
  }
  return(list(enrichment=Enrich,percentages=Percs,medianPercs=MedianPercs,maxPercs=maxPercs,minPercs=minPercs, pvalues=pvals))
}



rotate <- function(x) t(apply(x, 2, rev))


#######################
iterations <- 10000 # number of times the file was shuffled

stages<-c('Early Embryo', 'Larval Stage 3', 'Young Adult')
abrvs<-c('EE','L3','YA')

comparNames <- c('Heterochromatin'
                 ,'Repressed'
                 ,'Promoter'
                 ,'Transcribed'
                 ,"Repressed Enhancer"
                 ,"Active Enhancer"
                 ,'Low Signal')
reorderedInds<-c(7,1,2,5,6,3,4)

everything<-list()
atacEnrichs <- pvals <- atacPercs <- medianPercs <- minPercs <-  maxPercs<- matrix(nrow=length(comparNames), ncol=length(stages))

for (i in 1: length(stages)){
  stageAbrv <- abrvs[i]
  setwd(paste('/Users/acd13/Desktop/ATAC/Analysis/Enrichments/chromHMMState/enrichmentsWithNull/',stageAbrv,'/',stageAbrv,'_consensus_gDNAMasked',sep="")) # get to the right place
  
  # get the files
  atacEnrich<-list.files(pattern="Features.txt$") # these are ATAC-seq peaks near the features
  recips<-list.files(pattern="Input.txt$") # these are features near the ATAC-seq peaks
  
  #recips
  recips <- parseResults(recips[reorderedInds],comparNames[reorderedInds], stageAbrv)
  # ATAC
  main <- parseResults(atacEnrich[reorderedInds],comparNames[reorderedInds], stageAbrv)
  
  both <- list(reciprocal=recips, atac=main)  
  everything[[i]] <- both
  
  # because I rotate these they start in the columns, this is different than the gene based plotting
  atacEnrichs[,i] <- main$enrichment
  pvals [,i] <- main$pvalues
  atacPercs[,i] <- main$percentages
  medianPercs[,i] <- main$medianPercs
  minPercs[,i] <- main$minPercs
  maxPercs[,i] <- main$maxPercs
  
}


combined_enrich.rotated<-rotate(log2(atacEnrichs))
color=c('darkorchid4','goldenrod2','darkgreen')

# Plot the enrichment
setwd('/Users/acd13/Desktop/ATAC/Analysis/Enrichments/chromHMMState/enrichmentsWithNull')
pdf("enrichmentInStagePredictedStates.pdf", height=5, width=5)
par(mar=c(4.1,4.1,4.1,2.1))
barx <- barplot(combined_enrich.rotated[c(3,2,1),seq(7,1,-1)], 
                col=color[c(3,2,1)],
                xlab='Log2(enrichment over mappable genome)', 
                beside=T,
                yaxt='n',
                xlim=range(combined_enrich.rotated),
                horiz=T, cex.axis=1.2, cex.lab=1.2,
)

# the chromatin state nameswill need to be added separately
text(-1,26.5,"Transcribed")
text(-1,22.5,"Promoter")
text(-1,18.5,"Enhancer")
text(-1,14.5,"Repressed Enhnacer")
text(-1,10.5,"H3K27me3 Repressed")
text(1,6.5,"Heterochromatin")
text(1,2.5,"Low")
legend('topright',stages, col=color, pch=15, bty='n')
dev.off()



### Everything below here wasn't used

# Plot the percentages

combined_atacPercs.rotated<-rotate(atacPercs[c(-7,-8),])
combined_medianPercs.rotated<-rotate(medianPercs[c(-7,-8),])
library(scales)
color=c('darkorchid4',alpha('darkorchid4', 0.5),'goldenrod2',alpha('goldenrod2',0.5),'darkgreen',alpha('darkgreen',0.5))

enhancersOnlyToPlot <- rbind(combined_atacPercs.rotated[,4], combined_medianPercs.rotated[,4])
# Plot the enrichment
par(family='serif', cex=1.2)
par(mar=c(4.1,4.1,4.1,2.1))
barx <- barplot(rbind(atacPercs[3,],medianPercs[3,]), 
                col=color,
                beside=T
                ,ylim=c(0,35)
                ,ylab='% ATAC-seq peaks in enhancer predicted chromatin state'
                , cex.axis=1.2, cex.lab=1.2,
)
errorBarsTops <- c(atacPercs[3,1], maxPercs[3,1], atacPercs[3,2], maxPercs[3,2], atacPercs[3,3], maxPercs[3,3])
errorBarsBots <- c(atacPercs[3,1], minPercs[3,1], atacPercs[3,2], minPercs[3,2], atacPercs[3,3], minPercs[3,3])
# this genreates a lot of warnings, but no big deal
arrows(barx, errorBarsTops, barx, errorBarsBots, length=0.05, angle=90, code=3)

legend('top',stages, col=color[c(1,3,5)], pch=15, bty='n', cex=1.2)


