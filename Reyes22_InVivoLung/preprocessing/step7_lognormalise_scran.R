library(dplyr)
library(SingleCellExperiment)
library(zellkonverter)
library(scran)
library(scuttle)

reyes.sce <- readH5AD("../data/Reyes2022/reyes22_old_unspliced.h5ad", X_name="X")

assay(reyes.sce, "counts") <- assay(reyes.sce, "mature") + assay(reyes.sce, "ambiguous")
reyes.sce <- computePooledFactors(reyes.sce, clusters=reyes.sce$predictions_scanvi)
summary(sizeFactors(reyes.sce))

reyes.sce <- logNormCounts(reyes.sce)

assay(reyes.sce, "X") <- assay(reyes.sce, "logcounts")
assay(reyes.sce, "counts") <- NULL
assay(reyes.sce, "logcounts") <- NULL

writeH5AD(reyes.sce, "../data/Reyes2022/reyes22_old_unspliced.h5ad", X_name="X", compression="gzip")

