#!/usr/bin/python
"""This module runs data preprocessing."""

import numpy as np
import scanpy as sc
import seaborn as sns
from scipy.stats import median_abs_deviation

import json
import math
import multiprocessing
import os
import sys

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER = MOUNT_PATHS["input"]

# Setup of scanpy.
sc.settings.verbosity = 2
sc.settings.set_figure_params(
    dpi=80,
    facecolor="white",
    frameon=False,
)
sc.settings.figdir = MOUNT_PATHS["output"]

def is_outlier(data, metric: str, n_mad: int):
    """
    Detects outliers according to the defined metric by
    checking their distance from the median
    of the metric over the whole dataset.

    Keyword arguments:
    data -- the dataset
    metric -- the metric to use for outlier detection
    n_mad -- a datapoint is defined as outlier if its distance from 
    the metric median is more than n_mad times the 
    median absolute deviation (MAD)
    """
    m = data.obs[metric]
    return np.abs(m - np.median(m)) > n_mad * median_abs_deviation(m)

def process_data(file_path):
    """
    TODO
    """
    print(f"Processing file {file_path}")
    print("\tReading data...")
    adata = sc.read_10x_h5(filename = file_path)
    print("\tMaking variable names unique...")
    adata.var_names_make_unique()

    print("\tCalculating QC metrics...")
    # Marking mitochondrial genes.
    adata.var["mt"] = adata.var_names.str.lower().str.startswith("mt-")
    # Marking ribosomal genes.
    adata.var["ribo"] = adata.var_names.str.lower().str.startswith(("rps", "rpl"))
    # Marking hemoglobin genes.
    adata.var["hb"] = adata.var_names.str.lower().str.contains(("^hb[^(p)]"))

    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True
    )

    print("\tGenerating QC plots...")
    sns.displot(
        data = adata.obs["total_counts"],
        bins=100,
        kde=False
    ).set(
        xlabel ="Total counts per cell",
        ylabel = "Number of cells"
    ).savefig(f"{MOUNT_PATHS['output']}/histo_raw_total_counts.svg")
    sc.pl.violin(
        adata = adata,
        xlabel = "",
        ylabel = "Percentage of mitochondrial\ngene counts",
        keys="pct_counts_mt",
        show = False,
        save="_raw_mitochondiral_counts.svg"
    )
    n_genes_plot_raw = sc.pl.scatter(
        adata = adata,
        x = "total_counts",
        y = "n_genes_by_counts",
        color = "pct_counts_mt",
        title = "Percentage of mitochondrial gene counts",
        show = False
    )
    n_genes_plot_raw.set(
        xlabel ="Number of genes per cell",
        ylabel = "Total counts per cell"
    )
    n_genes_plot_raw.figure.savefig(f"{MOUNT_PATHS['output']}/scatter_raw_number_of_genes.svg")

    print("\tCalculating outliers...")
    adata.obs["outlier"] = (
        is_outlier(adata, "log1p_total_counts", 5)
        | is_outlier(adata, "log1p_n_genes_by_counts", 5)
        | is_outlier(adata, "pct_counts_in_top_20_genes", 5)
    )
    n_outliers = adata.obs.outlier[adata.obs.outlier].count()

    adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", 3) | (
        adata.obs["pct_counts_mt"] > 8
    )
    n_outliers_mt = adata.obs.mt_outlier[adata.obs.mt_outlier].count()

    print(f"\tNumber of cells before filtering: {adata.n_obs}")
    adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()

    print(f"\tNumber of cells after filtering: {adata.n_obs}")

    print("\tPlotting filtered data...")
    n_genes_plot_filtered = sc.pl.scatter(
        adata = adata,
        x = "total_counts",
        y = "n_genes_by_counts",
        color = "pct_counts_mt",
        title = "Percentage of mitochondrial gene counts",
        show = False
    )
    n_genes_plot_filtered.set(
        xlabel ="Number of genes per cell",
        ylabel = "Total counts per cell"
    )
    n_genes_plot_filtered.figure.savefig(f"{MOUNT_PATHS['output']}/scatter_filtered_number_of_genes.svg")

process_data(f"{INPUT_FOLDER}/sample_filtered_feature_bc_matrix.h5")
