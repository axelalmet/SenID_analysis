library(zellkonverter)
library(SingleCellExperiment)
library(DropletUtils)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(edgeR)
library(scuttle)

data_directory = '../data/Wechter2023/'

srr_runs = read.table(paste0(data_directory, 'SRR_Acc_List_static.txt'))$V1

knee_plot <- function(bc_rank) {
  knee_plt <- tibble(rank = bc_rank[["rank"]],
                     total = bc_rank[["total"]]) %>% 
    distinct() %>% 
    dplyr::filter(total > 0)
  annot <- tibble(inflection = metadata(bc_rank)[["inflection"]],
                  rank_cutoff = max(bc_rank$rank[bc_rank$total > metadata(bc_rank)[["inflection"]]]))
  p <- ggplot(knee_plt, aes(total, rank)) +
    geom_line() +
    geom_hline(aes(yintercept = rank_cutoff), data = annot, linetype = 2) +
    geom_vline(aes(xintercept = inflection), data = annot, linetype = 2) +
    scale_x_log10() +
    scale_y_log10() +
    annotation_logticks() +
    labs(y = "Rank", x = "Total UMIs")
  return(p)
}


empty_drops_plot <- function(empty_out, fdr) {
  empty_drops_plt <- tibble(total = empty_out$Total, 
                            logprob = -empty_out$LogProb,
                            is_cell = empty_out$FDR <= fdr_threshold)
  p <- ggplot(empty_drops_plt, aes(total, logprob, colour=is_cell)) +
    geom_point() +
    scale_x_log10() +
    scale_y_log10() +
    annotation_logticks() +
    labs(y = "-log(Probability)", x = "Total UMI count")
  return(p)
}

fdr_threshold <- 0.01

for (i in 1:length(srr_runs)) {
  run <- srr_runs[i]
  
  # Read in the file now
  file_directory <- paste0(data_directory, "out_", run, "_v3/counts_unfiltered/")
  
  sce <- readH5AD(paste0(file_directory, "adata_filtered.h5ad"))
  
  counts <- assay(sce, "mature") + assay(sce, "ambiguous") # Use spliced counts for analysis
  
  bc_rank <- barcodeRanks(counts)

  options(repr.plot.width=5, repr.plot.height=4)
  knee_plot(bc_rank)
  ggsave(paste0(file_directory, "knee_plot_", run, ".pdf"), width = 4, height = 4)
  
  set.seed(100*i)
    e.out <- emptyDrops(counts)
  is.cell <- (e.out$FDR <= fdr_threshold)
  print(table(Limited=e.out$Limited, Significant=is.cell))
  
  sum(is.cell, na.rm=TRUE)
  
  print(sum(is.cell, na.rm=TRUE))
  
  options(repr.plot.width=5, repr.plot.height=4)
  empty_drops_plot(e.out, fdr_threshold)
  ggsave(paste0(file_directory, "emptyDroplets_", run, ".pdf"), width = 4, height = 4)
  
  filtered.counts <- counts[,which(is.cell)]
  cells <- calculateAverage(filtered.counts)
  ambient.cells <- which(is.na(e.out$FDR))
  ambient <- rowSums(counts[, ambient.cells])
  
  options(repr.plot.width=5, repr.plot.height=4)
  res <- maPlot(ambient, cells, normalize = TRUE)
  dev.copy(pdf, paste0(file_directory, "MAplot_", run, ".pdf"), width = 4, height = 4)
  dev.off()

  sce <- sce[, which(is.cell)]

  print(dim(sce))
  
  writeH5AD(sce, paste0(file_directory, "adata_filtered.h5ad"), compression="gzip")
  
}

