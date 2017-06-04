args <- commandArgs(trailingOnly = TRUE)
heatmapRdataFile <- args[1]
allDirs <- strsplit(heatmapRdataFile, "/",fixed=T)
toName <- allDirs[[1]][(length(allDirs[[1]])-1)]
load(heatmapRdataFile)
g.order <- readRDS(args[2])


#set params
flood.q=.02 # inherited from original ngs plots

title2plot = toName
out.hm = paste(toName,'_cutomOrder.pdf',sep="")


RankNormalizeMatrix <- function(mat, low.cutoff) {
  # Rank-based normalization for a data matrix.
  # Args:
  #   mat: data matrix.
  #   low.cutoff: low value cutoff.
  # Return: rank normalized matrix.
  
  stopifnot(is.matrix(mat))
  
  concat.dat <- c(mat)
  low.mask <- concat.dat < low.cutoff
  concat.r <- rank(concat.dat)
  concat.r[low.mask] <- 0
  
  matrix(concat.r, nrow=nrow(mat))
}


SetupHeatmapDevice <- function(reg.list, uniq.reg, ng.list, pts, font.size=12,
                               unit.width=4, reduce.ratio=30) {
  # Configure parameters for heatmap output device. The output is used by 
  # external procedures to setup pdf device ready for heatmap plotting.
  # Args:
  #   reg.list: region list as in config file.
  #   uniq.reg: unique region list.
  #   ng.list: number of genes per heatmap in the order as config file.
  #   pts: data points (number of columns of heatmaps).
  #   font.size: font size.
  #   unit.width: image width per heatmap.
  #   reduce.ratio: how compressed are genes in comparison to data points? This 
  #                 controls image height.
  
  # Number of plots per region.
  reg.np <- sapply(uniq.reg, function(r) sum(reg.list==r))
  
  # Number of genes per region.
  reg.ng <- sapply(uniq.reg, function(r) {
    ri <- which(reg.list==r)[1]
    ng.list[ri]
  })
  
  # Adjustment ratio.
  origin.fs <- 12  # default font size.
  fs.adj.ratio <- font.size / origin.fs
  # Margin size (in lines) adjusted by ratio.
  m.bot <- fs.adj.ratio * 2
  m.lef <- fs.adj.ratio * 1.5
  m.top <- fs.adj.ratio * 2
  m.rig <- fs.adj.ratio * 1.5 
  key.in <- fs.adj.ratio * 1.0  # colorkey in inches.
  m.lef.diff <- (fs.adj.ratio - 1) * 1.5
  m.rig.diff <- (fs.adj.ratio - 1) * 1.5
  # Setup image size.
  hm.width <- (unit.width + m.lef.diff + m.rig.diff) * max(reg.np)
  ipl <- .2 # inches per line. Obtained from par->'mai', 'mar'.
  # Convert #gene to image height.
  reg.hei <- sapply(reg.ng, function(r) {
    c(key.in,  # colorkey + margin.
      r * unit.width / pts / reduce.ratio + 
        m.bot * ipl + m.top * ipl)  # heatmap + margin.
  })
  reg.hei <- c(reg.hei)
  hm.height <- sum(reg.hei)
  
  # Setup layout of the heatmaps.
  lay.mat <- matrix(0, ncol=max(reg.np), nrow=length(reg.np) * 2)
  fig.n <- 1  # figure plotting number.
  for(i in 1:length(reg.np)) {
    row.upper <- i * 2 - 1
    row.lower <- i * 2
    for(j in 1:reg.np[i]) {
      lay.mat[row.upper, j] <- fig.n;
      fig.n <- fig.n + 1
      lay.mat[row.lower, j] <- fig.n;
      fig.n <- fig.n + 1
    }
  }
  
  list(reg.hei=reg.hei, hm.width=hm.width, hm.height=hm.height, 
       lay.mat=lay.mat, heatmap.mar=c(m.bot, m.lef, m.top, m.rig) * ipl)
}


hd <- SetupHeatmapDevice(reg.list, uniq.reg, ng.list, pts, font.size, 
                         unit.width, rr)
reg.hei <- hd$reg.hei  # list of image heights for unique regions.
hm.width <- hd$hm.width  # image width.
hm.height <- hd$hm.height # image height.
lay.mat <- hd$lay.mat  # matrix for heatmap layout.
heatmap.mar <- hd$heatmap.mar # heatmap margins in inches.


# Setup basic parameters.
ncolor <- 256
if(bam.pair) {
  two.colors <- unlist(strsplit(hm.color, ':'))
  enrich.palette <- colorRampPalette(c('darkblue', 'white', 
                                     'red'), 
                                     bias=1, interpolate='spline')
} else {
  if(hm.color != "default") {
    enrich.palette <- colorRampPalette(c('snow', hm.color))
  } else {
    enrich.palette <- colorRampPalette(c('snow', 'red'))    
  }
}

hm_cols <- ncol(enrichList[[1]])

# Adjust X-axis tick position. In a heatmap, X-axis is [0, 1].
# Assume xticks$pos is from 0 to N(>0).
xticks$pos <- xticks$pos / tail(xticks$pos, n=1)  # scale to the same size.

# Define a function to calculate color breaks.
ColorBreaks <- function(max.e, min.e, bam.pair, ncolor) {
  # Args:
  #   max.e: maximum enrichment value to be mapped to color.
  #   min.e: minimum enrichment value to be mapped to color.
  #   bam.pair: boolean tag for bam-pair.
  #   ncolor: number of colors to use.
  # Returns: vector of color breaks.
  
  # If bam-pair is used, create breaks for positives and negatives 
  # separately. If log2 ratios are all positive or negative, use only 
  # half of the color space.
  if(bam.pair) {
    max.e <- ifelse(max.e > 0, max.e, 1)
    min.e <- ifelse(min.e < 0, min.e, -1)
    c(seq(min.e, 0, length.out=ncolor / 2 + 1),
      seq(0, max.e, length.out=ncolor / 2 + 1)[-1])
  } else {
    seq(min.e, max.e, length.out=ncolor + 1)
  }
}

if(grepl(",", color.scale)) {
  scale.pair <- unlist(strsplit(color.scale, ","))
  scale.min <- as.numeric(scale.pair[1])
  scale.max <- as.numeric(scale.pair[2])
  if(scale.min >= scale.max) {
    warning("Color scale min value is >= max value.\n")
  }
  flood.pts <- c(scale.min, scale.max)
  brk.use <- ColorBreaks(scale.max, scale.min, bam.pair, ncolor)
}

# If color scale is global, calculate breaks and quantile here.
if(color.scale == 'global') {
  flood.pts <- quantile(c(enrichList, recursive=T), c(flood.q, 1-flood.q))
  brk.use <- ColorBreaks(flood.pts[2], flood.pts[1], bam.pair, ncolor)
}

# Go through each unique region. 
# Do NOT use "dopar" in the "foreach" loops here because this will disturb
# the image order.
go.list <- vector('list', length(uniq.reg))
names(go.list) <- uniq.reg


pdf(out.hm, width=hm.width, height=hm.height, pointsize=font.size)
par(mai=heatmap.mar)
layout(lay.mat, heights=reg.hei)    

for(ur in uniq.reg) {
  # ur <- uniq.reg[1]
  plist <- which(reg.list==ur)  # get indices in the config file.
  
  # Combine all profiles into one.
  # enrichCombined <- do.call('cbind', enrichList[plist])
  enrichSelected <- enrichList[plist]
  
  # If color scale is region, calculate breaks and quantile here.
  if(color.scale == 'region') {
    flood.pts <- quantile(c(enrichSelected, recursive=T), 
                          c(flood.q, 1-flood.q))
    brk.use <- ColorBreaks(flood.pts[2], flood.pts[1], bam.pair, ncolor)
  }
  
  # Order genes.
  if(is.matrix(enrichSelected[[1]]) && nrow(enrichSelected[[1]]) > 1) {
    lowCutoffs <- v.low.cutoff[plist]
    
    # Order genes with a list of heatmap data.
    # Args: 
    #   enrichList: heatmap data in a list.
    #   lowCutoffs: low count cutoff for normalized count data.
    #   method: algorithm used to order genes.
    #   go.paras: gene ordering parameters.
    # Returns: a vector of gene orders. 
    
    rankList <- mapply(RankNormalizeMatrix,mat=enrichList, low.cutoff=lowCutoffs, SIMPLIFY=F)
    np <- length(enrichList)
    
    rankCombined <- do.call('cbind', rankList)
        
    go.list[[ur]] <- rev(rownames(enrichSelected[[1]][g.order, ]))
  } else {
    g.order <- NULL
    go.list[[ur]] <- g.order
  }
  
  
  # Go through each sample and do plot.
  for(pj in plist) {
    if(!is.null(g.order)) {
      enrichList[[pj]] <- enrichList[[pj]][g.order, ]
    }
    
    # If color scale is local, calculate breaks and quantiles here.
    if(color.scale == 'local') {
      flood.pts <- quantile(c(enrichList[[pj]], recursive=T), 
                            c(flood.q, 1-flood.q))
      brk.use <- ColorBreaks(flood.pts[2], flood.pts[1], bam.pair, 
                             ncolor)
    }
    
    # Flooding extreme values.
    enrichList[[pj]][ enrichList[[pj]] < flood.pts[1] ] <- flood.pts[1]
    enrichList[[pj]][ enrichList[[pj]] > flood.pts[2] ] <- flood.pts[2]
    
    # Draw colorkey.
    image(z=matrix(brk.use, ncol=1), col=enrich.palette(ncolor), 
          breaks=brk.use, axes=F, useRaster=T, main='Colorkey')
    axis(1, at=seq(0, 1, length.out=5), 
         labels=format(brk.use[seq(1, ncolor + 1, length.out=5)], 
                       digits=1), 
         lwd=1, lwd.ticks=1)
    
    # Draw heatmap.
    image(z=t(enrichList[[pj]]), col=enrich.palette(ncolor), 
          breaks=brk.use, axes=F, useRaster=T, main=title2plot[pj])
    
    axis(1, at=xticks$pos, labels=xticks$lab, lwd=1, lwd.ticks=1)
    
    
  }
}

out.dev <- dev.off()


