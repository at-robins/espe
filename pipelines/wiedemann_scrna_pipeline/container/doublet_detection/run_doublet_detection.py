#!/usr/bin/python
"""This module detects doublets."""

import anndata
import anndata2ri
import csv
import json
import logging
import numpy as np
import os
import scanpy as sc
import seaborn as sns
import rpy2.rinterface_lib.callbacks as rcb
import rpy2.robjects as ro
import warnings

from rpy2.robjects.packages import importr

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER = MOUNT_PATHS["dependencies"]["ambient_rna_removal"] + "/"

# Setup of rpy2.
rcb.logger.setLevel(logging.INFO)
ro.pandas2ri.activate()
anndata2ri.activate()
ro.r(
    '.libPaths(c("/usr/local/lib/R/site-library", "/usr/lib/R/site-library", "/usr/lib/R/library"))'
)


# Setup of scanpy.
sc.settings.verbosity = 2
sc.settings.set_figure_params(
    dpi=80,
    facecolor="white",
    frameon=False,
)
sc.settings.figdir = MOUNT_PATHS["output"]


def process_data(file_path_filtered, output_folder_path, metrics_writer):
    """
    Detects and marks doublets.
    """
    print(f"Processing file {file_path_filtered}", flush=True)
    print("\tReading filtered data...")
    adata_filtered = anndata.read_h5ad(file_path_filtered)

    print("\tLoading R dependencies...")
    importr("Seurat")
    importr("scater")
    importr("scDblFinder")
    # importr("BiocParallel")
    print("\tRunning scDblFinder...")
    doublet_detect_function = ro.r(
        """
        function(data) {
            set.seed(42)
            sce = scDblFinder(
                SingleCellExperiment(
                    list(counts=data),
                ) 
            )
            doublet_score = sce$scDblFinder.score
            doublet_class = sce$scDblFinder.class
            print(sce)
            print(doublet_score)
            print(doublet_class)
            return(sce)
        }
        """
    )
    doublet_detect_output = doublet_detect_function(
        adata_filtered.X.T,
    )
    print(doublet_detect_output)
    #     adata_filtered.layers["counts"] = adata_filtered.X
    #     adata_filtered.layers["soupX_counts"] = soupx_output.T
    #     adata_filtered.X = adata_filtered.layers["soupX_counts"]
    # n_cells_before_filter = adata_filtered.n_vars
    # print(f"Total number of features before filtering: {n_cells_before_filter}")
    # sc.pp.filter_genes(adata_filtered, min_cells=20)
    # n_cells_after_filter = adata_filtered.n_vars
    # print(f"Number of features after filtering: {n_cells_after_filter}")
    # print("\tWriting metrics to file...")
    # metrics_writer.writerow(
    #     [file_path_filtered, n_cells_before_filter, n_cells_after_filter]
    # )
    # print("\tWriting filtered data to file...")
    # adata_filtered.write(f"{output_folder_path}/corrected.h5ad", compression="gzip")

with open(
    f"{MOUNT_PATHS['output']}/metrics.csv", mode="w", newline="", encoding="utf-8"
) as csvfile:
    metrics_writer = csv.writer(
        csvfile, dialect="unix", delimiter=",", quotechar='"', quoting=csv.QUOTE_MINIMAL
    )
    metrics_writer.writerow(
        [
            "Sample",
            "Number of features before filtering",
            "Number of features after filtering",
        ]
    )
    # Iterates over all sample directories and processes them conserving the directory structure.
    for root, dirs, files in os.walk(INPUT_FOLDER):
        for file in files:
            if file.casefold().endswith("corrected.h5ad"):
                file_path_filtered = os.path.join(root, file)
                output_folder_path = os.path.join(
                    MOUNT_PATHS["output"], root.removeprefix(INPUT_FOLDER)
                )
                os.makedirs(output_folder_path, exist_ok=True)
                process_data(
                    file_path_filtered,
                    output_folder_path,
                    metrics_writer
                )
