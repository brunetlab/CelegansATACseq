library(pheatmap)

dds <- c('/Users/acd13/Desktop/ATAC/Analysis/Enrichments/TF_peaks_in_chromHMM/EE_bedtools_FisherTest_summaries/summits'
         ,'/Users/acd13/Desktop/ATAC/Analysis/Enrichments/TF_peaks_in_chromHMM/EE_bedtools_FisherTest_summaries'
         ,'/Users/acd13/Desktop/ATAC/Analysis/Enrichments/TFPeaks_in_other_peaks/eeTFPeaks'
         )
barplot_pdfs <- c('/Users/acd13/Desktop/ATAC/Analysis/Enrichments/TF_peaks_in_chromHMM/fishers_test_odds_ratios_TF100bpSummits_in_EEchromHMM.pdf'
                  ,'/Users/acd13/Desktop/ATAC/Analysis/Enrichments/TF_peaks_in_chromHMM/fishers_test_odds_ratios_TFPeaks_in_EEchromHMM.pdf'
                  ,'/Users/acd13/Desktop/ATAC/Analysis/Enrichments/TFPeaks_in_other_peaks/fishers_test_odds_ratios_TFPeaks_in_EETFPeaks.pdf'
                  )
heatmap_pdfs <- c('/Users/acd13/Desktop/ATAC/Analysis/Enrichments/TF_peaks_in_chromHMM/fishers_test_odds_ratios_TFPSummits_in_EEchromHMM_heatmaps.pdf',
                  '/Users/acd13/Desktop/ATAC/Analysis/Enrichments/TF_peaks_in_chromHMM/fishers_test_odds_ratios_TFPeaks_in_EEchromHMM_heatmaps.pdf',
                  '/Users/acd13/Desktop/ATAC/Analysis/Enrichments/TFPeaks_in_other_peaks/fishers_test_odds_ratios_TFPeaks_in_EETFPeaks_heatmaps.pdf'
                  )
name_inds <- c(11, 11, 10)
for (i in 1:length(dds)){
  data_dir <- dds[i]
  all_files <- list.files(path = data_dir, pattern="*_genicFisherEnrichSummary.txt_forPlotting$")
  data_list <- list()
  names <- list()
  for (f in all_files){
    if (grepl('DCC3160', f)) next
    name_parts <- strsplit(f, '_')
    strain <- sapply(name_parts, '[[', 1)
    gene <- sapply(name_parts, '[[', 2)
    ind <- length(data_list) + 1
    data_list[[ind]] <- read.table(paste(data_dir, f, sep="/"), header=F )
    names[[ind]] <- paste(strain, gene, sep="_")
  }
  
  odds_ratios <- sapply(data_list, '[[', 5)
  odds_ratios[is.infinite(odds_ratios)] <- NA
  rownames(odds_ratios) <- sapply(strsplit(sapply(data_list, '[[', 1)[,1], "/", fixed=T), '[[', name_inds[i])
  colnames(odds_ratios) <- unlist(names)
  colors <- rep(c('grey', 'wheat3'), ncol(odds_ratios)/2.0)
  colors[grep("EOR", colnames(odds_ratios))] <- 'green'
  colors[grep("ELT-3", colnames(odds_ratios))] <- 'firebrick'
  colors[grep("PHA-4", colnames(odds_ratios))] <- 'blue'
  colors[grep("HPL2", colnames(odds_ratios))] <- 'black'
  
  pdf(barplot_pdfs[i])
  for (j in 1:nrow(odds_ratios)){
    barx <- barplot(log2(odds_ratios[j,] + 1e-2)
            , main = rownames(odds_ratios)[j]
            , col=colors
            , xaxt='n'
            , ylab='Log2(Odds Ratio + 1e-2)'
            )
    
    axis(1
         , at=barx
         , labels = FALSE
         )
    text(barx-2
         , par("usr")[3] - 0.25
         , labels = colnames(odds_ratios)
         , srt = 45
         , pos = 1
         , xpd = TRUE
         , col=colors
         )
  }
  dev.off()
  pdf(heatmap_pdfs[i], width=12)
    pheatmap(log2(odds_ratios+ 1e-2))
  dev.off()
}
