#!/usr/bin/python
"""This module removes ambient RNA."""

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

from rpy2.robjects.packages import importr

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER = MOUNT_PATHS["dependencies"]["preprocessing"] + "/"

# Setup of rpy2.
rcb.logger.setLevel(logging.INFO)
ro.pandas2ri.activate()
anndata2ri.activate()
ro.r(".libPaths(c(\"/usr/local/lib/R/site-library\", \"/usr/lib/R/site-library\", \"/usr/lib/R/library\"))")


# Setup of scanpy.
sc.settings.verbosity = 2
sc.settings.set_figure_params(
    dpi=80,
    facecolor="white",
    frameon=False,
)
sc.settings.figdir = MOUNT_PATHS["output"]

def get_raw_file(raw_file_folder):
    """
    Returns the raw data file in the specified folder or None if
    none is present.

    Keyword arguments:
    raw_file_folder -- the directory path to search the raw file in
    """
    for file_path in os.listdir(raw_file_folder):
        full_path = os.path.join(raw_file_folder, file_path)
        if os.path.isfile(full_path) and file_path.casefold().endswith("raw_feature_bc_matrix.h5"):
            return full_path
    return None

def process_data(file_path_filtered, file_path_raw, output_folder_path, metrics_writer):
    """
    Removes ambient RNA.
    """
    print(f"Processing files {file_path_filtered} and {file_path_raw}")
    print("\tReading filtered data...")
    adata_filtered = anndata.read_h5ad(file_path_filtered)
    adata_raw = sc.read_10x_h5(file_path_raw)

    print("\tMaking variable names unique...")
    adata_raw.var_names_make_unique()

    print("\tNormalising data...")
    adata_filtered_tmp = adata_filtered.copy()
    sc.pp.normalize_per_cell(adata_filtered_tmp)
    sc.pp.log1p(adata_filtered_tmp)

    print("\tClustering...")
    sc.pp.pca(adata_filtered_tmp)
    sc.pp.neighbors(adata_filtered_tmp)
    sc.tl.leiden(adata_filtered_tmp, key_added="soupx_groups")
    soupx_groups = adata_filtered_tmp.obs["soupx_groups"]
    # Deletes the data reference to save memory.
    del adata_filtered_tmp
    soupx_cells = adata_filtered.obs_names
    soupx_genes = adata_filtered.var_names
    soupx_data_filtered = adata_filtered.X.T

    print("\tReading raw data...")
    adata_raw = sc.read_10x_h5(file_path_raw)
    soupx_data_raw = adata_raw.X.T
    # Deletes the data reference to save memory.
    del adata_raw

    print("\tLoading SoupX...")
    importr("SoupX")
    print("\tRunning SoupX...")
    soupx_function = ro.r('''
        function(data, data_tod, genes, cells, soupx_groups) {
            # specify row and column names of data
            rownames(data) = genes
            colnames(data) = cells
            # ensure correct sparse format for table of counts and table of droplets
            data <- as(data, "sparseMatrix")
            data_tod <- as(data_tod, "sparseMatrix")

            print(nrow(data))
            print(nrow(data_tod))
            # Generate SoupChannel Object for SoupX 
            sc = SoupChannel(data_tod, data, calcSoupProfile = FALSE)

            # Add extra meta data to the SoupChannel object
            soupProf = data.frame(row.names = rownames(data), est = rowSums(data)/sum(data), counts = rowSums(data))
            sc = setSoupProfile(sc, soupProf)
            # Set cluster information in SoupChannel
            sc = setClusters(sc, soupx_groups)

            # Estimate contamination fraction
            sc  = autoEstCont(sc, doPlot=FALSE)
            # Infer corrected table of counts and rount to integer
            out = adjustCounts(sc, roundToInt = TRUE)
        }
        ''')
    soupx_output = soupx_function(soupx_data_filtered, soupx_data_raw, soupx_genes, soupx_cells, soupx_groups)
    print(soupx_output)




with open(f"{MOUNT_PATHS['output']}/metrics.csv", mode="w", newline="") as csvfile:
    metrics_writer = csv.writer(
        csvfile, dialect="unix", delimiter=",", quotechar='"', quoting=csv.QUOTE_MINIMAL
    )
    metrics_writer.writerow(
        [
            "Sample",
            "Number of cells before filtering",
            "Number of cells after filtering",
            "Number of QC outlier cells",
            "Number of mitochondiral count outlier cells",
        ]
    )
    # Iterates over all sample directories and processes them conserving the directory structure.
    for root, dirs, files in os.walk(INPUT_FOLDER):
        for file in files:
            if file.casefold().endswith("preprocessed.h5ad"):
                file_path_filtered = os.path.join(root, file)
                folder_path_raw = os.path.join(
                    MOUNT_PATHS["input"], root.removeprefix(INPUT_FOLDER)
                )
                output_folder_path = os.path.join(
                    MOUNT_PATHS["output"], root.removeprefix(INPUT_FOLDER)
                )
                os.makedirs(output_folder_path, exist_ok=True)
                process_data(file_path_filtered, get_raw_file(folder_path_raw), output_folder_path, metrics_writer)
