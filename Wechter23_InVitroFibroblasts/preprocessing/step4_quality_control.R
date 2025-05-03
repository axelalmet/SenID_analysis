library(zellkonverter)
library(SingleCellExperiment)
library(DropletUtils)
library(scuttle)

data_directory = '../data/Wechter2023/'

srr_runs = read.table(paste(data_directory, 'SRR_Acc_List_static.txt', sep=''))$V1

for (i in 1:length(srr_runs)) {
  run <- srr_runs[i]
  
  # Read in the file now
  file_directory <- paste(data_directory, "out_", run, "_v3/counts_unfiltered/", sep="")
  
  sce <- readH5AD(paste(file_directory, "adata_filtered.h5ad", sep=""))

  assay(sce, "counts") <- assay(sce, "mature") + assay(sce, "ambiguous")
  
  is.mito <- grep("^MT-", rownames(sce))
  
  per.cell <- perCellQCMetrics(sce, subsets=list(Mito=is.mito))
  
  print(summary(per.cell$sum))
  is.outlier <- isOutlier(per.cell$sum, type="both", nmads=3, log=TRUE)
  print(attr(is.outlier, "threshold"))
  
  is.mito.outlier <- isOutlier(per.cell$subsets_Mito_percent, type="higher", nmads=3, log=TRUE)
  print(attr(is.mito.outlier, "threshold"))

  sce <- sce[, (!is.outlier)&(!is.mito.outlier)]
  print(dim(sce))
  
  writeH5AD(sce, paste(file_directory, "adata_filtered.h5ad", sep=""), compression="gzip")
  
}
  
