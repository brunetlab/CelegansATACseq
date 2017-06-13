### Note that only the TF peaks in other peaks was reported

library(pheatmap)

# set up input and output names
dds <- c('/Users/acd13/Desktop/ATAC/Analysis/Enrichments/TF_peaks_in_chromHMM/L3_bedtools_FisherTest_summaries/summits'
         ,'/Users/acd13/Desktop/ATAC/Analysis/Enrichments/TFPeaks_in_other_peaks/l3TFPeaks'
         )
heatmap_pdfs <- c('/Users/acd13/Desktop/ATAC/Analysis/Enrichments/TF_peaks_in_chromHMM/fishers_test_odds_ratios_TFPSummits_in_L3chromHMM_heatmaps.pdf',
                  '/Users/acd13/Desktop/ATAC/Analysis/Enrichments/TFPeaks_in_other_peaks/fishers_test_odds_ratios_TFPeaks_in_L3TFPeaks_heatmaps.pdf'
                  )

name_ind <- 11

for (i in 1:length(dds)){
  data_dir <- dds[i]
  all_files <- list.files(path = data_dir, pattern="*_genicFisherEnrichSummary.txt_forPlotting$")
  data_list <- list()
  names <- list()
  for (f in all_files){
    # DCC3160 is not a trusted data source
    if (grepl('DCC3160', f)) next
    # set up the name
    name_parts <- strsplit(f, '_')
    strain <- sapply(name_parts, '[[', 1)
    gene <- sapply(name_parts, '[[', 2)
    ind <- length(data_list) + 1
    data_list[[ind]] <- read.table(paste(data_dir, f, sep="/"), header=F )
    names[[ind]] <- paste(strain, gene, sep="_")
  }
  # pull out the relevant data
  odds_ratios <- sapply(data_list, '[[', 5)
  # Set the diagnol (self comparisons) to NA
  odds_ratios[is.infinite(odds_ratios)] <- NA
  # update names
  rownames(odds_ratios) <- sapply(strsplit(sapply(data_list, '[[', 1)[,1], "/", fixed=T), '[[', name_ind)
  colnames(odds_ratios) <- unlist(names)
  
  # a small value was added to avoid log2 of 0
  pdf(heatmap_pdfs[i], width=12)
    pheatmap(log2(odds_ratios+ 1e-2))
  dev.off()
}
