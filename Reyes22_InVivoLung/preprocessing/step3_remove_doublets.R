library(zellkonverter)
library(SingleCellExperiment)
library(DropletUtils)
library(scDblFinder)

data_directory = '../data/Reyes2022/'

srr_runs_gfppos = read.table(paste0(data_directory, 'SRR_Acc_List_Old_GFPPos.txt'))$V1
srr_runs_gfpneg = read.table(paste0(data_directory, 'SRR_Acc_List_Old_GFPNeg.txt'))$V1

srr_runs <- c(srr_runs_gfppos, srr_runs_gfpneg)

for (i in 1:length(srr_runs)) {
  run <- srr_runs[i]
  
  # Read in the file now
  file_directory <- paste0(data_directory, "Old_GFPPos/out_", run, "_v3/counts_unfiltered/")
  if (run %in% srr_runs_gfpneg)
  {
    file_directory <- paste0(data_directory, "Old_GFPNeg/out_", run, "_v3/counts_unfiltered/")
  }
  
  sce <- readH5AD(paste0(file_directory, "adata_filtered.h5ad"))
  assay(sce, "counts") <- assay(sce, "mature") + assay(sce, "ambiguous") + assay(sce, "nascent")
  
  print(dim(sce))
  sce <- scDblFinder(sce)
  sce <- sce[, sce$scDblFinder.class != 'doublet']
  print(dim(sce))
  
  writeH5AD(sce, paste0(file_directory, "adata_filtered.h5ad"), compression="gzip")
  
}

