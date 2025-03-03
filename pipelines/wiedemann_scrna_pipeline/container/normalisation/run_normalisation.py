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
import scipy as sp
import rpy2.rinterface_lib.callbacks as rcb
import rpy2.robjects as ro
import warnings

from matplotlib import pyplot as plt
from rpy2.robjects.packages import importr
from scipy.sparse import csr_matrix
from scipy.sparse import issparse

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER = next(iter(MOUNT_PATHS["dependencies"].values()))

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


def shifted_logarithm(data):
    """
    Normalises the data using the shifted logarithm algorithm.
    Useful for dimensionality reduction and identification of
    differentialy present features.
    """
    print("\tApplying shifted logarithm normalisation...")
    scales_counts = sc.pp.normalize_total(data, target_sum=None, inplace=False)
    data.layers["log1p_norm"] = sc.pp.log1p(scales_counts["X"], copy=True)


def analytic_pearson_residuals(data):
    """
    Normalises the data using the analytic Pearson residuals algorithm.
    Useful for identification of variable features and rare cell populations.
    """
    print("\tApplying analytic Pearson residuals normalisation...")
    analytic_pearson = sc.experimental.pp.normalize_pearson_residuals(
        data, inplace=False
    )
    data.layers["analytic_pearson_residuals"] = csr_matrix(analytic_pearson["X"])


def proportional_fitting(data_matrix):
    """
    Performs proportional fitting on the input matrix.
    """
    sum_array = data_matrix.sum(axis=1).A.ravel()
    array_mean = sum_array.mean()
    return sp.sparse.diags(array_mean / sum_array) @ data_matrix


def pf_log1p_pf(data):
    """
    Normalises the data using the proportional-fitting-log1p-proportional-fitting algorithm.
    Useful if the sequencing depth is very different between samples.
    """
    print("\tApplying PFlog1pPF normalisation...")
    data.layers["PFlog1pPF"] = csr_matrix(
        proportional_fitting(np.log1p(proportional_fitting(data.X)))
    )


def scran(data):
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
    sns.histplot(data.obs["n_counts"], bins=100, kde=False, ax=axes_1)
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
    fig.savefig(os.path.join(output_folder_path, "histo_normalised_total_counts.svg"))
    plt.close(fig)


input_file_path = os.path.join(INPUT_FOLDER, "merged.h5ad")
print(f"Processing file {input_file_path}", flush=True)
print("\tReading data...")
output_folder_path = MOUNT_PATHS["output"]
adata_merged = anndata.read_h5ad(input_file_path)
shifted_logarithm(adata_merged)
analytic_pearson_residuals(adata_merged)
scran(adata_merged)
pf_log1p_pf(adata_merged)
plot_normalised_data(adata_merged, output_folder_path)
print("\tWriting normalised data to file...")
os.makedirs(output_folder_path, exist_ok=True)
adata_merged.write(
    os.path.join(output_folder_path, "normalised.h5ad"), compression="gzip"
)
