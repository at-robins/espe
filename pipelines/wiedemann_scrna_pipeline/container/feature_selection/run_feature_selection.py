#!/usr/bin/python
"""This module selects relevant features."""

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

from matplotlib import pyplot as plt
from rpy2.robjects.packages import importr
from scipy.sparse import csr_matrix
from scipy.sparse import issparse

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER = MOUNT_PATHS["dependencies"]["normalisation"] + "/"
BATCHED_SUBFOLDER = "batched"
UNBATCHED_SUBFOLDER = "unbatched"
INPUT_FOLDER_BATCHED = f"{INPUT_FOLDER}{BATCHED_SUBFOLDER}/"
INPUT_FOLDER_UNBATCHED = f"{INPUT_FOLDER}{UNBATCHED_SUBFOLDER}/"
TOP_FEATURES = 4000
DEFAULT_BATCH_KEY = "replicate_name"

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
    scanpy=True,
    # # In case of bitmap exports, use high quality.
    dpi=300,
    dpi_save=300,
    # Export as SVG.
    format="svg",
    vector_friendly=False,
    # Use transparent background.
    transparent=True,
    facecolor=None,
    # Remove frames.
    frameon=False,
)
sc.settings.figdir = MOUNT_PATHS["output"]


def is_batch(adata: anndata.AnnData) -> bool:
    """
    Returns true if the data contains batch information.
    """
    return (
        DEFAULT_BATCH_KEY in adata.obs_keys()
        and adata.obs[DEFAULT_BATCH_KEY].cat.categories.size > 1
    )


def calculate_deviance_unbatched(adata: anndata.AnnData):
    """
    Returns binomial deviance.
    """
    print("\tLoading R dependencies...")
    importr("scry")
    print("\tRunning R code...")
    scry_function = ro.r(
        """
        function(data) {
            sce <- devianceFeatureSelection(data, assay="X")
            return(
                rowData(sce)$binomial_deviance
            )
        }
        """
    )
    return scry_function(adata).T


def calculate_deviance_batched(adata: anndata.AnnData):
    """
    Returns binomial deviance for batched data.
    """
    batches = adata.obs[DEFAULT_BATCH_KEY].cat.categories
    sum_squared_deviance = np.zeros(adata.n_vars)
    for batch in batches:
        print(f"\tCalculating deviance for batch {batch}...")
        adata_subset = adata[adata.obs[DEFAULT_BATCH_KEY] == batch]
        # Calculate deviance and replace features not present in the 
        # subset with 0.
        deviance = np.nan_to_num(calculate_deviance_unbatched(adata_subset))
        # Sum up the square root of the batch deviances to negatively 
        # impact features not present in all batches.
        sum_squared_deviance = sum_squared_deviance + np.sqrt(deviance)
    return sum_squared_deviance


def process_data(file_path_input, output_folder_path):
    """
    Select relevant features.
    """
    print(f"Processing file {file_path_input}", flush=True)
    print("\tReading filtered data...")
    adata = anndata.read_h5ad(file_path_input)

    print("\tCalculating highly deviant genes...")
    if is_batch(adata):
        print("\tFound batch information.", flush=True)
        batch_key = DEFAULT_BATCH_KEY
        binomial_deviance = calculate_deviance_batched(adata)
    else:
        print("\tNo batch information found. Proceeding without...", flush=True)
        batch_key = None
        binomial_deviance = calculate_deviance_unbatched(adata)

    sorted_indices = binomial_deviance.argsort()
    if sorted_indices.size > TOP_FEATURES:
        print(f"\t{sorted_indices.size} features present. Using top {TOP_FEATURES}.")
        sorted_indices = sorted_indices[-TOP_FEATURES:]
    else:
        print(f"\tOnly {sorted_indices.size} features present. Using all of them.")
    mask = np.zeros(adata.var_names.shape, dtype=bool)
    mask[sorted_indices] = True
    adata.var["highly_deviant"] = mask
    adata.var["binomial_deviance"] = binomial_deviance
    print("\tCalculating highly variable genes...")
    sc.pp.highly_variable_genes(adata, layer="scran_normalisation", batch_key=batch_key)

    print("\tPlotting data...")
    fig, ax = plt.subplots()
    sns.scatterplot(
        data=adata.var, x="means", y="dispersions", hue="highly_deviant", s=5, ax=ax
    )
    ax.set_xlim(None, 4)
    ax.set_ylim(None, 4)
    ax.set(xlabel="Means", ylabel="Dispersion")
    legend = ax.get_legend()
    legend.get_texts()[0].set_text("No")
    legend.get_texts()[1].set_text("Yes")
    legend.set_title("Highly deviant")
    fig.tight_layout()
    fig.savefig(f"{output_folder_path}/scatter_dispersion.svg")

    print("\tWriting filtered data to file...")
    adata.write(
        f"{output_folder_path}/filtered_feature_bc_matrix.h5ad", compression="gzip"
    )


# Iterates over all batched sample directories and processes them conserving the directory structure.
for root, dirs, files in os.walk(INPUT_FOLDER_BATCHED):
    for file in files:
        if file.casefold().endswith("filtered_feature_bc_matrix.h5ad"):
            file_path_input = os.path.join(root, file)
            output_folder_path = os.path.join(
                MOUNT_PATHS["output"], root.removeprefix(INPUT_FOLDER)
            )
            os.makedirs(output_folder_path, exist_ok=True)
            process_data(file_path_input, output_folder_path)

# Iterates over all unbatched sample directories and processes them conserving the directory structure.
for root, dirs, files in os.walk(INPUT_FOLDER_UNBATCHED):
    for file in files:
        if file.casefold().endswith("filtered_feature_bc_matrix.h5ad"):
            file_path_input = os.path.join(root, file)
            output_folder_path = os.path.join(
                MOUNT_PATHS["output"], root.removeprefix(INPUT_FOLDER)
            )
            os.makedirs(output_folder_path, exist_ok=True)
            process_data(file_path_input, output_folder_path)
