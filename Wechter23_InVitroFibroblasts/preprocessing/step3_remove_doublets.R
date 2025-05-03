library(zellkonverter)
library(SingleCellExperiment)
library(DropletUtils)
library(scDblFinder)

data_directory = '../data/Wechter2023/'

srr_runs = read.table(paste(data_directory, 'SRR_Acc_List_static.txt', sep=''))$V1

for (i in 1:length(srr_runs)) {
  run <- srr_runs[i]
  
  # Read in the file now
  file_directory <- paste(data_directory, "out_", run, "_v3/counts_unfiltered/", sep="")
  
  sce <- readH5AD(paste(file_directory, "adata_filtered.h5ad", sep=""))
  assay(sce, "counts") <- assay(sce, "mature") + assay(sce, "ambiguous")
  
  print(dim(sce))
  
  sce <- scDblFinder(sce)
  sce <- sce[, sce$scDblFinder.class != 'doublet']
  
  print(dim(sce))
  writeH5AD(sce, paste(file_directory, "adata_filtered.h5ad", sep=""), compression="gzip")
  
}