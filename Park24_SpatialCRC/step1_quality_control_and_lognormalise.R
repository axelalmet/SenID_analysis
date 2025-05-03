library(SpotSweeper)
library(SpaNorm)
library(zellkonverter)
library(SingleCellExperiment)
library(SpatialExperiment)
library(scuttle)
library(SpaNorm)

# Save the data using this solution proposed here: https://github.com/theislab/zellkonverter/issues/56
spe_to_h5ad <- function(spe, file){
  stopifnot(all(rownames(colData(spe)) == rownames(spatialCoords(spe))))
  colData(spe) <- cbind(colData(spe), spatialCoords(spe))
  writeH5AD(spe, file, X_name="logcounts", compression="gzip")
}

data_directory <- "../data/Park2024/"

sample_ids <- c("P1", "P2", "P3", "P4")

for (sample_id in sample_ids)

{
    samples <- file.path(data_directory, sample_id, "outs")
    names(samples) <- sample_id

    list.files(samples[1])
    list.files(file.path(samples[1], "spatial"))
    file.path(samples[1], "filtered_feature_bc_matrix")

    file.path(file.path(samples[1], "spatial"), "scalefactors_json.json")

    spe <- read10xVisium(samples=samples,
                        type = "HDF5", data="filtered",
                        images = "hires", load = TRUE)
    colnames(spe) <- make.unique(colnames(spe)) # Making the column names unique.
    rownames(spe) <- make.unique(rowData(spe)$symbol)


    # Sanity checkt o make sure we remove spots not in tissue
    spe <- spe[, spe$in_tissue == 1]

    # identifying the mitochondrial transcripts
    is.mito <- rownames(spe)[grepl("^MT-", rownames(spe))]

    sessionInf# Add a counts assay (right now it's just the X assay)

    # calculating QC metrics for each spot using scuttle
    spe <- addPerCellQCMetrics(spe, subsets = list(Mito = is.mito))
    colnames(colData(spe))

    # library size
    spe <- localOutliers(spe,
                        metric = "sum",
                        direction = "lower",
                        log = TRUE
    )

    # unique genes
    spe <- localOutliers(spe,
                        metric = "detected",
                        direction = "lower",
                        log = TRUE
    )

    # mitochondrial percent
    spe <- localOutliers(spe,
                        metric = "subsets_Mito_percent",
                        direction = "higher",
                        log = FALSE
    )

    # combine all outliers into "local_outliers" column
    spe$local_outliers <- as.logical(spe$sum_outliers) |
    as.logical(spe$detected_outliers) |
    as.logical(spe$subsets_Mito_percent_outliers)

    # Remove outliers
    spe <- spe[, spe$local_outliers == FALSE]
    spe <- spe[, colSums(counts(spe)) != 0]

    # Variance stabilise
    spe = SpaNorm(spe)

    spe_to_h5ad(spe, file=paste0(dir, "/park24_", sample_id, "_spe.h5ad"))

}
