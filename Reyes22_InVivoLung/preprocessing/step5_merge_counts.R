library(zellkonverter)
library(SingleCellExperiment)

data_directory = '../data/Reyes2022/'

srr_runs_gfppos = read.table(paste0(data_directory, 'SRR_Acc_List_Old_GFPPos.txt'))$V1
srr_runs_gfpneg = read.table(paste0(data_directory, 'SRR_Acc_List_Old_GFPNeg.txt'))$V1

srr_runs <- c(srr_runs_gfppos, srr_runs_gfpneg)

sce.list <- list()
genes.list <- list()

for (i in 1:length(srr_runs)) {
  run <- srr_runs[i]
  
  # Read in the file now
  file_directory <- paste0(data_directory, "Old_GFPPos/out_", run, "_v3/counts_unfiltered/")
  if (run %in% srr_runs_gfpneg)
  {
    file_directory <- paste0(data_directory, "Old_GFPNeg/out_", run, "_v3/counts_unfiltered/")
  }
  
  sce <- readH5AD(paste0(file_directory, "adata_filtered.h5ad"))
  
  sce.list[[i]] <- sce
  genes.list[[i]] <- rownames(sce)

}

genes.all <- Reduce(intersect, genes.list)
print(length(genes.all))

for (i in 1:length(sce.list))
{
  sce.list[[i]] <- sce.list[[i]][genes.all, ]
  rowData(sce.list[[i]])[, 'scDblFinder.selected'] <- NULL
  colData(sce.list[[i]])[, 'n_cells'] <- NULL
  rowData(sce.list[[i]])[, 'n_cells'] <- NULL
}

sce.merged <- Reduce(SingleCellExperiment::cbind, sce.list)

print(dim(sce.merged))
assay(sce, "counts") <- NULL

writeH5AD(sce.merged, X_name="X", paste0(data_directory, "reyes22_old_unspliced.h5ad"), compression="gzip")
  
