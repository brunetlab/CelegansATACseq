stages <- c('EE', 'L3', 'YA')
d <- list()
for (stage in stages){
  f <- paste0('/Users/acd13/Desktop/ATAC/Analysis/Enrichments/',stage,'ChromHMMStateAnnotatedPeaks/metaPeaks_',stage,'ChromHMMStateAnnotationSummary.txt')
  d[[length(d) +1]] <- read.table(f, header=F, sep="\t")
}

all <- data.frame(EE = d[[1]][-1,2], L3 = d[[2]][-1,2], YA=d[[3]][-1,2]) # get rid of the total
rownames(all) <- c('Low signal', 'Active enhancer', 'Repressed', 'Heterochromatin', 'TSS', 'Gene body', 'Repressed enhancer')
reorder_ind <- c(1,4,3,7,2,5,6)
colors <- c('gray88','midnightblue', 'darkgrey', 'darkkhaki', 'darkorange2', 'firebrick3', 'green3')
pdf('~/Dropbox/stacked_chromHMM_counts_consensus_atacPeaks.pdf', width=5, height=5)
barplot(as.matrix(all[reorder_ind,])
        , ylim=c(0,35000)
        , ylab = 'Number of ATAC-seq consensus peaks'
        , las = 1
        , col = colors
        )
dev.off()

### To compare to the above I also calcualted the total number of 100bp bins that each state in each stage covered
# this was done using this line:
#for stage in EE L3 YA; do while read line ; do echo -e ${line} $(cat $line | awk '{sum += ($3 - $2)/100} END {print sum}'); done < /Users/acd13/Desktop/ATAC/Analysis/Enrichments/chromHMMState/enrichmentsWithNull/${stage}/${stage}MergedStates_finalSelect.txt ; done > ~/Desktop/ATAC/Analysis/checkingChromHMMs/bins_totals.tsv
dfile <- '~/Desktop/ATAC/Analysis/checkingChromHMMs/bins_totals.tsv'
dat <- read.table(dfile)
d <- list()

for (i in 1:length(stages)){
  stage <- stages[i]
  stage_dat <- vector(mode='numeric', length=length(reorder_ind))
  for (j in 1:length(reorder_ind)){
    stage_dat[j] <- dat[(i-1)*length(reorder_ind) + j,2]
  }
  d[[length(d) +1]] <- stage_dat
}

all <- data.frame(EE = d[[1]], L3 = d[[2]], YA=d[[3]]) # get rid of the total
rownames(all) <- c('Low signal', 'Active enhancer', 'Repressed', 'Heterochromatin', 'TSS', 'Gene body', 'Repressed enhancer')
pdf('~/Dropbox/stacked_chromHMM_bin_counts.pdf', width=5, height=5)
barplot(as.matrix(all[reorder_ind,])
        , ylab = 'Number of 100bp bins'
        , las = 1
        , col = colors
)
dev.off()
