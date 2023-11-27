#!/usr/bin/python
"""This module performs batch correction."""

import anndata
import anndata2ri
import csv
import json
import logging
import math
import multiprocessing
import numpy as np
import os
import scanpy as sc
import seaborn as sns
import rpy2.rinterface_lib.callbacks as rcb
import rpy2.robjects as ro
import warnings

from rpy2.robjects.packages import importr

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER = MOUNT_PATHS["dependencies"]["feature_selection"] + "/"
BATCHED_SUBFOLDER = "batched"
INPUT_FOLDER_BATCHED = f"{INPUT_FOLDER}{BATCHED_SUBFOLDER}/"

DEFAULT_BATCH_KEY = "batch"

# Setup of rpy2.
rcb.logger.setLevel(logging.INFO)
ro.pandas2ri.activate()
anndata2ri.activate()
ro.r(
    '.libPaths(c("/usr/local/lib/R/site-library", "/usr/lib/R/site-library", "/usr/lib/R/library"))'
)


# Setup of scanpy.
sc.settings.verbosity = 2

threads = math.floor(multiprocessing.cpu_count() * 0.8)
if threads < 1:
    threads = 1


def process_data(input_file_path, output_folder_path):
    """
    Corrects batch effects.
    """
    print(f"Processing file {input_file_path}", flush=True)
    print("\tReading data...")
    # TODO: REMOVE!
    if not "tumour" in input_file_path:
        return
    adata = anndata.read_h5ad(input_file_path)
    adata = adata[:, adata.var["highly_deviant"]].copy()
    print(f"\tUsing data {adata.n_vars} features for {adata.n_obs} cells")
    sciber_genes = adata.var_names
    sciber_cells = []
    sciber_data = []
    batches = adata.obs[DEFAULT_BATCH_KEY].cat.categories

    for batch in batches:
        print(f"\tExtracting batch {batch}...")
        adata_subset = adata[adata.obs[DEFAULT_BATCH_KEY] == batch]
        sciber_data.append(adata_subset.X.T)
        sciber_cells.append(adata_subset.obs_names)

    print("\tLoading R dependencies...")
    importr("Seurat")
    importr("SCIBER")
    print("\tRunning batch correction...")
    sciber_function = ro.r(
        """
            function(data, genes, cells, threads) {
                # Constructs dataframes.
                for(i in 1:length(data)){
                    rownames(data[[i]]) = genes
                    colnames(data[[i]]) = cells[[i]]
                    data[[i]] <- as.matrix(as(data[[i]], "sparseMatrix"))
                }
                return(SCIBER(input_batches = data, n_core = threads, seed = 42))
            }
            """
    )
    sciber_output = sciber_function(sciber_data, sciber_genes, sciber_cells, threads)
    print(sciber_output.keys())
    print(type(sciber_output))

    # print("\tUpdating data with doublet information...")
    # doublet_detect_output = doublet_detect_function(
    #     adata_filtered.X.T,
    # )
    # doublet_classes = doublet_detect_output.obs["scDblFinder.class"].to_numpy()
    # adata_filtered.obs["doublet_score"] = doublet_detect_output.obs[
    #     "scDblFinder.score"
    # ].to_numpy()
    # adata_filtered.obs["doublet_class"] = doublet_classes

    # print("\tWriting metrics to file...")
    # metrics_writer.writerow(
    #     [
    #         input_file_path.removeprefix(INPUT_FOLDER),
    #         adata_filtered.n_obs,
    #         np.count_nonzero(doublet_classes == "singlet"),
    #         np.count_nonzero(doublet_classes == "doublet"),
    #     ]
    # )
    # print("\tWriting filtered data to file...")
    # adata.write(
    #     f"{output_folder_path}/filtered_feature_bc_matrix.h5ad", compression="gzip"
    # )


# Iterates over all sample directories and processes them conserving the directory structure.
for root, dirs, files in os.walk(INPUT_FOLDER_BATCHED):
    for file in files:
        if file.casefold().endswith("filtered_feature_bc_matrix.h5ad"):
            file_path_input = os.path.join(root, file)
            output_folder_path = os.path.join(
                MOUNT_PATHS["output"], root.removeprefix(INPUT_FOLDER_BATCHED)
            )
            os.makedirs(output_folder_path, exist_ok=True)
            process_data(file_path_input, output_folder_path)
