# originally written 20Jun2016 by Aaron Daugherty - Brunet Lab, Stanford University
# this program basically calclulates all pairwise jaccards for a set of genomic locations
# in this case gene based as well as reapts and the region of interest, gDNA peaks from N2 gDNA-ATAC (Input control)

# set up
Sys.setenv(PATH = "$PATH:/usr/bin:/usr/local/bin:/Users/acd13/Softwares/bedtools2-2.21.0/bin/")
library(pheatmap)
library(pvclust)
# get files
peaksDir <- '/Users/acd13/Desktop/ATAC/Analysis/gDNAPeaks/peaksToUse'
bedFiles <- list.files(peaksDir, pattern="*.bed$")
peakFiles <- list.files(peaksDir, pattern="*Peak$")
filesToUse <- c(bedFiles, peakFiles)

names <-c(
 "TTS1kbWindow"
# ,"ATAC" # I had initially included ATAC peaks, having forgotten that these peaks had been filtered for exactly these gDNA peaks
 ,"Repeats"
 ,"Distal"
 ,"Exons"
 ,"Extended Promoter"
 ,"Genes body"
 ,"Introns"
# ,"H3K9me3" # everything else was genomic location based, so I nixed this one
 ,"gDNA peaks"
)

jaccards.nucleotides <- matrix(0, nrow=length(filesToUse),ncol=length(filesToUse))
rownames(jaccards.nucleotides) <- colnames(jaccards.nucleotides) <- names
# run all of the pairwise comparisons by using system calls to Bedtools
setwd(peaksDir)
for (i in 1:length(filesToUse)){
  for (j in 1:length(filesToUse)){
    # bedtools does it counting the number of bp overlap
    rawResult <- system(paste("bedtools jaccard -a ", filesToUse[i], " -b ", filesToUse[j],sep=""), intern=TRUE)
    jaccards.nucleotides[i,j] <- as.numeric(strsplit(rawResult[2], "\t")[[1]][3])    
  }
}

# output the results
setwd("..")
#write.table(jaccards.nucleotides, file="gDNAJaccardVsGenomicLoci.txt", quote=F)
clustered.nuc<-pvclust(jaccards.nucleotides, method.hclust="complete",nboot=10000)

pdf("~/Dropbox/gDNANucleotideJaccardLog10.pdf", height=5, width=5)
  barplot(jaccards.nucleotides[,ncol(jaccards.nucleotides)]
          , ylab='Bedtools Jaccard Index')
  barplot(jaccards.nucleotides[-ncol(jaccards.nucleotides),ncol(jaccards.nucleotides)]
          , main='without self comparison'
          , ylab='Bedtools Jaccard Index')
  pheatmap(log10(jaccards.nucleotides + 1e-4))
  plot(clustered.nuc)
dev.off()

### I did not use the below, it was just for experimententing
#
## The above shows that these repeats are overlapping repeats the most.
# I'm intrigued to see if there is a particular type of repeat that it overlaps
# I'll test that here
#
peaksDir <- '/Users/acd13/Desktop/ATAC/Analysis/gDNAPeaks/peaksToUse/subtypesOfRepeats/'
filesToUse <- list.files(peaksDir, pattern="*.bed$")
# this adds all of the repeats subtypes, but I also want to add the gDNA peaks, and for copmarison, I'll include all of the repats together
filesToUse[length(filesToUse)+1] <- '/Users/acd13/Desktop/ATAC/Analysis/gDNAPeaks/peaksToUse/ce10_repeatMasker.bed'
filesToUse[length(filesToUse)+1] <- '/Users/acd13/Desktop/ATAC/Analysis/gDNAPeaks/peaksToUse/N2_gDNA.filt.nodup.adjusted.macs2.1q5e-2_peaks.narrowPeak'
# since this is taking a while, I want to go ahead and include everything I might be interested in.
# in this case that includes the ATAC-seq peaks
filesToUse[length(filesToUse)+1] <- '/Users/acd13/Desktop/ATAC/Analysis/gDNAPeaks/peaksToUse/allRepsAllButL1MetaPeaks_smartMerged300bpSummitDist.bedNotUsing'


jaccards.nucleotides <- matrix(0, nrow=length(filesToUse),ncol=length(filesToUse))
rownames(jaccards.nucleotides) <- colnames(jaccards.nucleotides) <- filesToUse
# run all of the pairwise comparisons by using system calls to Bedtools
setwd(peaksDir)
for (i in 1:length(filesToUse)){
  for (j in 1:length(filesToUse)){
    # bedtools does it counting the number of bp overlap
    rawResult <- system(paste("bedtools jaccard -a ", filesToUse[i], " -b ", filesToUse[j],sep=""), intern=TRUE)
    jaccards.nucleotides[i,j] <- as.numeric(strsplit(rawResult[2], "\t")[[1]][3])
  }
}

# output the results
setwd("../..")
saveRDS(jaccards.nucleotides, file="gDNAAndMetaATACSeqPeaksJaccardVsRepeatSubTypes.rds")
write.table(jaccards.nucleotides, file="gDNAAndMetaATACSeqPeaksJaccardVsRepeatSubTypes.txt", quote=F)
clustered.nuc<-pvclust(jaccards.nucleotides, method.hclust="complete",nboot=10000)

pdf("~/Dropbox/gDNAAndMetaATACSeqPeaksNucleotideJaccardLog10_repeatSubtypes.pdf", height=25, width=25) # I know this is huge, but I wanted to be able see things

  # first view a pheatmap and clustering of all vs all
  pheatmap(log10(jaccards.nucleotides + 1e-4, main="all repeats, gDNA and Meta ATAC peaks")
  plot(clustered.nuc, main='all repeats, gDNA and Meta ATAC peaks') 
  
  # then view the jaccard index of all repeats vs the ATAC meta peaks, with and without a self comparison
  colOfInterest <- ncol(jaccards.nucleotides)
  barplot(jaccards.nucleotides[,colOfInterest], ylab='Bedtools Jaccard Index', main="All repeats, and gDNA peaks vs meta peaks) 
  barplot(jaccards.nucleotides[-colOfInterest,colOfInterest], ylab='Bedtools Jaccard Index', main="All repeats, and gDNA peaks vs meta peaks, no self comparison) 
  barplot(jaccards.nucleotides[-colOfInterest,colOfInterest], ylab='Bedtools Jaccard Index', main="All repeats, and gDNA peaks vs meta peaks, no self comparison) 
  
  # do the same, but for the gDNA peaks
  colOfInterest <- ncol(jaccards.nucleotides)-1
  barplot(jaccards.nucleotides[,colOfInterest], ylab='Bedtools Jaccard Index', main="All repeats, and meta peaks vs gDNA peaks) 
  barplot(jaccards.nucleotides[-colOfInterest,colOfInterest], ylab='Bedtools Jaccard Index', main="All repeats, and meta peaks vs gDNA peaks, no self comparison) 
  #
  #
dev.off()

