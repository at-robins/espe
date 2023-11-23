#!/usr/bin/python
"""This module normalises data."""

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
INPUT_FOLDER_BATCHED = MOUNT_PATHS["dependencies"]["integration"] + "/"
INPUT_FOLDER_UNBATCHED = MOUNT_PATHS["dependencies"]["doublet_detection"] + "/"

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


def shifted_logarithm(data, output_folder_path):
    """
    Normalises the data using the shifted logarithm algorithm.
    Useful for dimensionality reduction and identification of
    differentialy present features.
    """
    print("\tApplying shifted logarithm normalisation...")
    scales_counts = sc.pp.normalize_total(data, target_sum=None, inplace=False)
    data.layers["log1p_norm"] = sc.pp.log1p(scales_counts["X"], copy=True)


def analytic_pearson_residuals(data, output_folder_path):
    """
    Normalises the data using the analytic Pearson residuals algorithm.
    Useful for identification of variable features and
    """
    print("\tApplying analytic Pearson residuals normalisation...")
    analytic_pearson = sc.experimental.pp.normalize_pearson_residuals(
        data, inplace=False
    )
    data.layers["analytic_pearson_residuals"] = csr_matrix(analytic_pearson["X"])


def scran(data, output_folder_path):
    """
    Normalises the data using the scran algorithm.
    Useful for batch correction.
    """
    print("\tApplying scran normalisation...")
    adata_pp = data.copy()
    sc.pp.normalize_total(adata_pp)
    sc.pp.log1p(adata_pp)
    sc.pp.pca(adata_pp, n_comps=15)
    sc.pp.neighbors(adata_pp)
    sc.tl.leiden(adata_pp, key_added="groups")
    data_mat = adata_pp.X.T
    input_groups = adata_pp.obs["groups"]
    # Converts the matrix if possible. See: https://github.com/MarioniLab/scran/issues/70
    if issparse(data_mat):
        if data_mat.nnz > 2**31 - 1:
            data_mat = data_mat.tocoo()
        else:
            data_mat = data_mat.tocsc()
    del adata_pp

    print("\tLoading R dependencies...")
    importr("scran")
    importr("BiocParallel")
    print("\tRunning R code...")
    scran_function = ro.r(
        """
        function(data_mat, input_groups) {
            return(
                sizeFactors(
                    computeSumFactors(
                        SingleCellExperiment(
                            list(counts=data_mat)
                        ), 
                        clusters = input_groups,
                        min.mean = 0.1,
                        BPPARAM = MulticoreParam()
                    )
                )
            )
        }
        """
    )
    size_factors = scran_function(data_mat, input_groups)
    print("\tUpdating data with scran information...")
    data.obs["size_factors"] = size_factors
    scran = data.X / data.obs["size_factors"].values[:, None]
    data.layers["scran_normalisation"] = csr_matrix(sc.pp.log1p(scran))


def plot_normalised_data(data, output_folder_path):
    """
    Plots the normalised data.
    """
    print("\tPlotting normalised data...")
    fig, ((axes_1, axes_2), (axes_3, axes_4)) = plt.subplots(2, 2, figsize=(10, 10))

    axes_1.set_title("Before normalisation")
    sns.histplot(data.obs["total_counts"], bins=100, kde=False, ax=axes_1)
    axes_1.set(xlabel="Total counts per cell", ylabel="Number of cells")

    axes_2.set_title("Shifted logarithm")
    sns.histplot(data.layers["log1p_norm"].sum(1), bins=100, kde=False, ax=axes_2)
    axes_2.set(xlabel="Total normalised counts per cell", ylabel="Number of cells")
    axes_2.get_legend().remove()

    axes_3.set_title("Analytic Pearson residuals")
    sns.histplot(
        data.layers["analytic_pearson_residuals"].sum(1),
        bins=100,
        kde=False,
        ax=axes_3,
    )
    axes_3.set(xlabel="Total normalised counts per cell", ylabel="Number of cells")
    axes_3.get_legend().remove()

    axes_4.set_title("Scran")
    sns.histplot(
        data.layers["scran_normalisation"].sum(1),
        bins=100,
        kde=False,
        ax=axes_4,
    )
    axes_4.set(xlabel="Total normalised counts per cell", ylabel="Number of cells")
    axes_4.get_legend().remove()

    fig.tight_layout()
    fig.savefig(f"{output_folder_path}/histo_normalised_total_counts.svg")


def normalise_data(file_path_filtered, output_folder_path):
    """
    Normalises the data using different algorithms.
    """
    print(f"Processing file {file_path_filtered}", flush=True)
    print("\tReading filtered data...")
    adata_filtered = anndata.read_h5ad(file_path_filtered)
    shifted_logarithm(adata_filtered, output_folder_path)
    analytic_pearson_residuals(adata_filtered, output_folder_path)
    scran(adata_filtered, output_folder_path)
    plot_normalised_data(adata_filtered, output_folder_path)
    print("\tWriting filtered data to file...")
    adata_filtered.write(
        f"{output_folder_path}/filtered_feature_bc_matrix.h5ad", compression="gzip"
    )


def process_data(input_directory, output_directory):
    """
    Iterates over all sample directories and processes the according files
    while conserving the directory structure.
    """
    for root, dirs, files in os.walk(input_directory):
        for file in files:
            if file.casefold().endswith("filtered_feature_bc_matrix.h5ad"):
                file_path_filtered = os.path.join(root, file)
                output_folder_path = os.path.join(
                    MOUNT_PATHS["output"],
                    output_directory,
                    root.removeprefix(input_directory),
                )
                os.makedirs(output_folder_path, exist_ok=True)
                normalise_data(file_path_filtered, output_folder_path)


process_data(INPUT_FOLDER_UNBATCHED, "unbatched")
process_data(INPUT_FOLDER_BATCHED, "batched")
