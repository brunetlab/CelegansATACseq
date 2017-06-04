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
 ,"Repeats"
 ,"Distal"
 ,"Exons"
 ,"Extended Promoter"
 ,"Genes body"
 ,"Introns"
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
dev.off()
