library('scales')
library(matrixStats)
setwd('~/Dropbox')
dds <- c('/Users/acd13/Desktop/ATAC/Analysis/tfs/qualityTfs/insertSizeStatsAtPeaks/ATACSignalSummitsInTFs/L3/rawResults/perTfPeak/insertSizes/L3/',
         "/Users/acd13/Desktop/ATAC/Analysis/tfs/qualityTfs/insertSizeStatsAtPeaks/TFPeakSummits/L3/rawResults/perTfPeak/insertSizes/L3/"
         )
histPdfs <- c('allL3TFPeaksInsertSizeHist_allpeak.pdf',
              'allL3TFPeaksInsertSizeHist_justSummits.pdf'
              )
aggregatehistPdfs <- c('allL3TFPeaksInsertSizeAggregateHist_allpeak.pdf',
              'allL3TFPeaksInsertSizeAggregateHist_justSummits.pdf'
                )
averagedhistPdfs <- c('allL3TFPeaksInsertSizeAvgHist_allpeak.pdf',
                       'allL3TFPeaksInsertSizeAvgHist_justSummits.pdf'
                )
ecdf_pdfs <- c('allL3TFPeaksInsertSizeECDFPlots_justSummits_onlyTFs.pdf',
               'allL3TFPeaksInsertSizeECDFPlots_wholePeaks_onlyTFs.pdf'
               )
all_ecdf_pdfs <- c('allL3TFPeaksInsertSizeECDFPlots_justSummits.pdf',
               'allL3TFPeaksInsertSizeECDFPlots_wholePeaks.pdf'
)
for (i in 1:length(dds)){
  data_dir <- dds[i] 
  
  toload <- list.files(data_dir, pattern = "*.txt$")
  allTF_hs <- list()
  nonTF_hs <- list()
  
  
  pdf(histPdfs[i], width =5, height =5)
  for (f in toload){
      if (grepl('DCC3160', f)) next
      d <- read.table(paste0(data_dir,f), sep="\t", header=F)[,2]
      h <- hist(d
       #, breaks=100
       , breaks=seq(0,2000,4)       
       , xlim = c(0,500)
        , xlab = 'Median insert size (bp)'
        , col = 'goldenrod2'
        , las =1
        , main = gsub(".txt","",f)
        , freq = FALSE
        )
      h$density <- h$counts/sum(h$counts)*100
      if (grepl('OP81_EOR',f)){
          eorh <- h
      }else if(grepl('SDQ2354_HDA1',f) | grepl('SDQ2340_HPL2',f)){
          nonTF_hs[[length(nonTF_hs) + 1]] <- h
      }else{
          allTF_hs[[length(allTF_hs) + 1]] <- h
      }
    }
  dev.off()
  pdf(aggregatehistPdfs[i], width=5, height=5)
  plot(nonTF_hs[[1]]
       , xlim=c(0,500)
       , ylim=c(0,11)
       , freq=F
       , col=alpha('blue',0.25)
       ,las=1
       ,ylab='Percentage of peaks'
       ,xlab='Median ATAC-seq insert size in TF peak'
       ,lty='blank'
  )
  plot(nonTF_hs[[2]]
       , freq=F
       , col=alpha('blue',0.2)
       ,lty='blank'
       , add=T)
  for (j in 1:length(allTF_hs)){
    plot(allTF_hs[[j]]
         , freq=F
         , col=alpha('black',0.1)
         ,lty='blank'
         , add=T)
  }  
  plot(eorh, freq=F, col=alpha('red',0.5),lty='blank', add=T)
  dev.off()  
    
  nonTF_density <- rbind(nonTF_hs[[1]]$density, nonTF_hs[[2]]$density)
  
  allTF_density <- rbind(allTF_hs[[1]]$density, allTF_hs[[2]]$density)
  
  for (j in 3:length(allTF_hs)){
    allTF_density <- rbind(allTF_density,allTF_hs[[j]]$density)
    
  }
  pdf(averagedhistPdfs[i], width=5, height=5)
  plot(allTF_hs[[1]]$mids
       ,colMeans(allTF_density)
       , xlim=c(0,500)
       , type ='h'
       , lwd = 3
       ,ylab='Percentage of peaks'
       ,xlab='Median ATAC-seq insert size in TF peaks (6bp bins)'
       ,col=alpha('black', 0.5)
       ,las=1
       )
  lines(allTF_hs[[1]]$mids
        ,colMeans(nonTF_density)
        , lwd = 3
        , type ='h'
       ,col=alpha('blue',0.5)
  )
  lines(eorh$mids
        ,eorh$density
        ,col=alpha('red',0.5)
        , lwd = 3
        , type ='h'
  )
  dev.off()
  
  eor1_ind <- grep("OP81_EOR-1.txt", toload)
  hpl2 <- grep("SDQ2340_HPL2.txt", toload)
  hda <- grep("SDQ2354_HDA1.txt", toload)
  dcc <- grep("DCC3160_EOR-1.txt", toload)
  others <- toload[-1 * c(hpl2,hda,eor1_ind,dcc)]
  
  pdf(ecdf_pdfs[i])
  plot(ecdf(read.table(paste0(data_dir,toload[eor1_ind])
                       , sep="\t", header=F)[,2])
       , col='white'
       , xlim=c(0,600)
       , main = 'Cumulative Distributions of ATAC-seq median insert sizes in L3 TFs'     
       , las = 1
       , ylab = 'Cumulative distibution portion'
       , xlab='Median ATAC-seq insert size (bp)'
       )
  
  for (f in others){
    lines(ecdf(read.table(paste0(data_dir,f), sep="\t", header=F)[,2]), col=alpha('black', 0.1))
  }
  lines(ecdf(read.table(paste0(data_dir,toload[eor1_ind])
                        , sep="\t", header=F)[,2]), col=alpha('red', 0.4))
  legend('right',,legend=c('Other TFs', 'EOR-1'), col=c('black', 'red'), pch=16)
  dev.off()
  
  noTfs <- toload[c(hpl2,hda)]
  pdf(all_ecdf_pdfs[i])
  plot(ecdf(read.table(paste0(data_dir,toload[eor1_ind])
                       , sep="\t", header=F)[,2])
       , col='white'
       , xlim=c(0,600)
       , main = 'Cumulative Distributions of ATAC-seq median insert sizes in L3 TFs'     
       , las = 1
       , ylab = 'Cumulative distibution portion'
       , xlab='Median ATAC-seq insert size (bp)'
  )
  
  for (f in others){
    lines(ecdf(read.table(paste0(data_dir,f), sep="\t", header=F)[,2]), col=alpha('black', 0.2))
  }
  for (f in noTfs){
    lines(ecdf(read.table(paste0(data_dir,f), sep="\t", header=F)[,2]), col=alpha('blue', 0.2))
  }
  lines(ecdf(read.table(paste0(data_dir,toload[eor1_ind])
                        , sep="\t", header=F)[,2]), col=alpha('red', 0.4))
  legend('right',,legend=c('Other TFs', 'HDA-1/HPL-2', 'EOR-1'), col=c('black', 'blue', 'red'), pch=16)
  dev.off()
}