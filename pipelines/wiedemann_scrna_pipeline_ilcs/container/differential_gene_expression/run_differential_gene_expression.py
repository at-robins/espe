#!/usr/bin/python
"""This module detects doublets."""

import anndata
import anndata2ri
import csv
import json
import logging
import numpy as np
import os
import pandas as pd
import scanpy as sc
import seaborn as sns
import rpy2.rinterface_lib.callbacks as rcb
import rpy2.robjects as ro
import warnings

from rpy2.robjects.packages import importr

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER = MOUNT_PATHS["dependencies"]["ilc_composition"] + "/"
CELL_TYPE_KEY = "cell_type"
BATCH_KEY = "batch"

# Setup of rpy2.
rcb.logger.setLevel(logging.INFO)
ro.pandas2ri.activate()
anndata2ri.activate()
ro.r(
    '.libPaths(c("/usr/local/lib/R/site-library", "/usr/lib/R/site-library", "/usr/lib/R/library"))'
)


# Setup of scanpy.
sc.settings.verbosity = 2


def process_data(file_path_filtered, output_folder_path, metrics_writer):
    """
    Detects and marks doublets.
    """
    print(f"Processing file {file_path_filtered}", flush=True)
    print("\tReading data...")
    adata = anndata.read_h5ad(file_path_filtered)

    cell_types = adata.obs[CELL_TYPE_KEY].cat.categories
    batches = adata.obs[BATCH_KEY].cat.categories
    pseudobulk_dataframe = pd.DataFrame()

    print(f"\tNumber of total observations: {adata.n_obs}")

    for cell_type in cell_types:
        for batch in batches:
            print(f"\tProcessing batch {batch} and cell type {cell_type}")
            cell_type_mask = adata.obs[CELL_TYPE_KEY] == cell_type
            batch_mask = adata.obs[BATCH_KEY] == batch
            adata_subset = adata[np.logical_and(cell_type_mask, batch_mask)]
            aggregated_row = adata_subset.to_df().agg(np.sum)
            aggregated_row["celltype"] = cell_type
            aggregated_row["replicate"] = batch
            pseudobulk_dataframe = pd.concat(
                [pseudobulk_dataframe, aggregated_row.to_frame().T],
                ignore_index=False,
                join="outer",
            )

    pseudobulk_dataframe.fillna(0, inplace=True)
    print(pseudobulk_dataframe)

    # print("\tLoading R dependencies...")
    # importr("Seurat")
    # importr("edgeR")
    # print("\tRunning R...")
    # doublet_detect_function = ro.r(
    #     """
    #     function(data){
    #         # create an edgeR object with counts and grouping factor
    #         y <- DGEList(assay(data, "X"), group = colData(data)$label)
    #         # filter out genes with low counts
    #         print("Dimensions before subsetting:")
    #         print(dim(y))
    #         print("")
    #         keep <- filterByExpr(y)
    #         y <- y[keep, , keep.lib.sizes=FALSE]
    #         print("Dimensions after subsetting:")
    #         print(dim(y))
    #         print("")
    #         # normalize
    #         y <- calcNormFactors(y)
    #         # create a vector that is concatentation of condition and cell type that we will later use with contrasts
    #         group <- paste0(colData(data)$label, ".", colData(data)$cell_type)
    #         replicate <- colData(data)$replicate
    #         # create a design matrix: here we have multiple donors so also consider that in the design matrix
    #         design <- model.matrix(~ 0 + group + replicate)
    #         # estimate dispersion
    #         y <- estimateDisp(y, design = design)
    #         # fit the model
    #         fit <- glmQLFit(y, design)
    #         return(list("fit"=fit, "design"=design, "y"=y))
    #     }
    #     """
    # )
    # print("\tUpdating data with doublet information...")
    # doublet_detect_output = doublet_detect_function(
    #     adata.X.T,
    # )
    # doublet_classes = doublet_detect_output.obs["scDblFinder.class"].to_numpy()
    # adata.obs["doublet_score"] = doublet_detect_output.obs[
    #     "scDblFinder.score"
    # ].to_numpy()
    # adata.obs["doublet_class"] = doublet_classes

    # print("\tWriting metrics to file...")
    # metrics_writer.writerow(
    #     [
    #         file_path_filtered.removeprefix(INPUT_FOLDER),
    #         adata.n_obs,
    #         np.count_nonzero(doublet_classes == "singlet"),
    #         np.count_nonzero(doublet_classes == "doublet"),
    #     ]
    # )
    # print("\tWriting filtered data to file...")
    # adata.write(
    #     f"{output_folder_path}/filtered_feature_bc_matrix.h5ad", compression="gzip"
    # )


with open(
    f"{MOUNT_PATHS['output']}/metrics.csv", mode="w", newline="", encoding="utf-8"
) as csvfile:
    metrics_writer = csv.writer(
        csvfile, dialect="unix", delimiter=",", quotechar='"', quoting=csv.QUOTE_MINIMAL
    )
    metrics_writer.writerow(
        [
            "Sample",
            "Total number of cells",
            "Number of singlets",
            "Number of doublets",
        ]
    )
    # Iterates over all sample directories and processes them conserving the directory structure.
    for root, dirs, files in os.walk(INPUT_FOLDER):
        for file in files:
            if file.casefold().endswith("filtered_feature_bc_matrix.h5ad"):
                file_path_filtered = os.path.join(root, file)
                output_folder_path = os.path.join(
                    MOUNT_PATHS["output"], root.removeprefix(INPUT_FOLDER)
                )
                os.makedirs(output_folder_path, exist_ok=True)
                process_data(file_path_filtered, output_folder_path, metrics_writer)
