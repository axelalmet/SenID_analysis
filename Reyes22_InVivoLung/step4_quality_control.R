library(zellkonverter)
library(SingleCellExperiment)
library(DropletUtils)
library(scuttle)


data_directory = '/Users/axelalmet/Documents/CellCellCommunicationModelling/Data/Reyes2022/'

srr_runs_gfppos = read.table(paste(data_directory, 'SRR_Acc_List_Old_GFPPos.txt', sep=''))$V1
srr_runs_gfpneg = read.table(paste(data_directory, 'SRR_Acc_List_Old_GFPNeg.txt', sep=''))$V1

srr_runs <- c(srr_runs_gfppos, srr_runs_gfpneg)

for (i in 1:length(srr_runs)) {
  run <- srr_runs[i]
  
  # Read in the file now
  file_directory <- paste(data_directory, "Old_GFPPos/out_", run, "_v3/counts_unfiltered/")
  if (run %in% srr_runs_gfpneg)
  {
    file_directory <- paste(data_directory, "Old_GFPNeg/out_", run, "_v3/counts_unfiltered/")
  }
  
  sce <- readH5AD(paste(file_directory, "adata_filtered.h5ad"))
  assay(sce, "counts") <- assay(sce, "mature") + assay(sce, "ambiguous")
  
  is.mito <- grep("^mt-", rownames(sce))
  
  per.cell <- perCellQCMetrics(sce, subsets=list(Mito=is.mito))
  
  print(summary(per.cell$sum))
  is.outlier <- isOutlier(per.cell$sum, type="both", nmads=3, log=TRUE)
  print(attr(is.outlier, "threshold"))
  
  is.mito.outlier <- isOutlier(per.cell$subsets_Mito_percent, type="higher", nmads=3, log=TRUE)
  print(attr(is.mito.outlier, "threshold"))
  
  print(dim(sce))
  sce <- sce[, (!is.outlier)&(!is.mito.outlier)]
  print(dim(sce))
  
  
  writeH5AD(sce, paste(file_directory, "adata_filtered.h5ad"), compression="gzip")
  
}


