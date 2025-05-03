library(dplyr)
library(SingleCellExperiment)
library(zellkonverter)
library(scran)
library(scuttle)

wechter.sce <- readH5AD("../data/Wechter2023/wechter23_perturbed_unspliced.h5ad", X_name="X")

# In case we deleted the "counts"
assay(wechter.sce, "counts") <- assay(wechter.sce, "mature") + assay(wechter.sce, "ambiguous") 

clusters <- quickCluster(wechter.sce)
wechter.sce <- computePooledFactors(wechter.sce, clusters=clusters)
summary(sizeFactors(wechter.sce))

wechter.sce <- logNormCounts(wechter.sce)
assay(wechter.sce, "counts") <- NULL
writeH5AD(wechter.sce, "../data/Wechter2023/wechter23_perturbed_unspliced.h5ad", X_name="X", compression="gzip")

