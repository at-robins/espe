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
TOP_FEATURES = 4000

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
    dpi=300,
    facecolor="white",
    frameon=False,
)
sc.settings.figdir = MOUNT_PATHS["output"]

def process_data(file_path_input, output_folder_path):
    """
    Select relevant features.
    """
    print(f"Processing file {file_path_input}", flush=True)
    print("\tReading filtered data...")
    adata = anndata.read_h5ad(file_path_input)

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
    binomial_deviance = scry_function(adata).T

    print("\tUpdating data...")
    idx = binomial_deviance.argsort()
    if idx.size > TOP_FEATURES:
        print(f"\t{idx.size} features present. Using top {TOP_FEATURES}.")
        idx = idx[-TOP_FEATURES:]
    else:
        print(f"\tOnly {idx.size} features present. Using all of them.")
    mask = np.zeros(adata.var_names.shape, dtype=bool)
    mask[idx] = True
    adata.var["highly_deviant"] = mask
    adata.var["binomial_deviance"] = binomial_deviance
    sc.pp.highly_variable_genes(adata, layer="scran_normalisation")

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
    adata.write(f"{output_folder_path}/feature_selection.h5ad", compression="gzip")


# Iterates over all sample directories and processes them conserving the directory structure.
for root, dirs, files in os.walk(INPUT_FOLDER):
    for file in files:
        if file.casefold().endswith("normalised.h5ad"):
            file_path_input = os.path.join(root, file)
            output_folder_path = os.path.join(
                MOUNT_PATHS["output"], root.removeprefix(INPUT_FOLDER)
            )
            os.makedirs(output_folder_path, exist_ok=True)
            process_data(file_path_input, output_folder_path)
