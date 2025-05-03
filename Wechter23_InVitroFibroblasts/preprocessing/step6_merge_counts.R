library(zellkonverter)
library(SingleCellExperiment)
library(SEtools)

data_directory = '/Users/axelalmet/Documents/CellCellCommunicationModelling/Data/Wechter2023/'

srr_runs = read.table(paste0(data_directory, 'SRR_Acc_List_static.txt'))$V1

sce.list <- list()
genes.list <- list()

for (i in 1:length(srr_runs)) {
  run <- srr_runs[i]

  # Read in the file now
  file_directory <- paste0(data_directory, "UnsplicedSpliced/out_", run, "_v3/counts_unfiltered/")

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
assay(sce.merged, "counts") <- NULL

print(dim(sce.merged))
writeH5AD(sce.merged, X_name="X", paste0(data_directory, "wechter23_merged_unspliced.h5ad"), compression="gzip")

