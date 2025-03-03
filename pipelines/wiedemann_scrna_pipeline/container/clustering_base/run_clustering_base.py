#!/usr/bin/python
"""This module performs data clustering."""

import anndata
import csv
import json
import logging
import numpy as np
import os
import pandas as pd
import pathvalidate
import rpy2.robjects as ro
import scanpy as sc
import seaborn as sns

from matplotlib import pyplot as plt
from rpy2.robjects.packages import importr

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER = next(iter(MOUNT_PATHS["dependencies"].values()))
MINIMUM_CELL_NUMBER = 20
PCS_MINIMUM_INFERENCE = 50
PCS_MINIMUM_NEIGHBOURS = 30


# The replicates.
REPLICATE_KEY = "replicate_name"
# The samples that consists of different replicates.
SAMPLE_TYPE_KEY = "sample_type"
# The clustering information.
CLUSTER_KEY = "leiden_clustering"
# The sample groups that have been clustered together.
CLUSTERING_GROUP_KEY = "clustering_group"

# Setup of scanpy.
sc.settings.verbosity = 2
sc.settings.set_figure_params(
    scanpy=True,
    # # In case of bitmap exports, use high quality.
    dpi=300,
    dpi_save=300,
    # Export as SVG.
    format="svg",
    vector_friendly=True,
    # Use transparent background.
    transparent=True,
    facecolor=None,
    # Remove frames.
    frameon=False,
)
sc.settings.figdir = ""


file_path_input = os.path.join(INPUT_FOLDER, "feature_selection.h5ad")
output_folder_path = MOUNT_PATHS["output"]
os.makedirs(output_folder_path, exist_ok=True)

print(f"Processing file {file_path_input}...")
print("\tReading data...", flush=True)
adata = anndata.read_h5ad(file_path_input)

adata.var["highly_variable"] = adata.var["highly_deviant"]

# Ensures that there are no missing values from the integration.
np.nan_to_num(adata.layers["log1p_norm"], copy=False, nan=0.0)
adata.X = adata.layers["log1p_norm"]

print(f"\tPerforming clustering for {adata.n_obs} cells...", flush=True)
if adata.n_obs < MINIMUM_CELL_NUMBER:
    print("\tNot enough cells. Skipping sample...")
else:
    sc.pp.pca(adata, svd_solver="arpack", use_highly_variable=True)
    n_pcs_max = adata.obsm["X_pca"].shape[1]
    if n_pcs_max < PCS_MINIMUM_INFERENCE:
        n_pcs_inference = n_pcs_max
    else:
        n_pcs_inference = PCS_MINIMUM_INFERENCE

    # Plots PCA variance.
    sc.pl.pca_variance_ratio(
        adata,
        n_pcs=n_pcs_inference,
        show=True,
        save=False,
    )
    plt.savefig(
        os.path.join(
            output_folder_path,
            "pca_variance.svg",
        )
    )
    plt.close()

    print("\tCalculating optimal number of PCs to use...", flush=True)
    pca_variances = adata.uns["pca"]["variance"]
    pca_sds = np.sqrt(pca_variances)
    # Sorts in decreasing order as required by findPC.
    pca_sds[::-1].sort()
    print("\tLoading R dependencies...")
    importr("findPC")
    print("\tRunning R...")
    find_pc_function = ro.r(
        """
        function(data, max_pca){
            return(findPC(sdev = data, number = max_pca))
        }
        """
    )
    # R returns a matrix with single float value so some conversion is needed.
    inferred_number_of_pcs = int(
        np.asmatrix(
            find_pc_function(
                ro.FloatVector(pca_sds.tolist()),
                n_pcs_inference,
            )
        )[0, 0]
    )
    print(f"\tInferred number of PCs: {inferred_number_of_pcs}", flush=True)

    # The inferred number of PCs might be too low,
    # so a minimum of PCs to use for clustering
    # was defined.
    if n_pcs_max < PCS_MINIMUM_NEIGHBOURS:
        pcs_neighbours = n_pcs_max
    else:
        pcs_neighbours = PCS_MINIMUM_NEIGHBOURS

    if inferred_number_of_pcs >= pcs_neighbours:
        pcs_neighbours = inferred_number_of_pcs
    else:
        print(
            f"\tInferred PC number is smaller than the defined minimum and is thus ignored."
        )

    print(f"\tUsing {pcs_neighbours} PCs for clustering...", flush=True)

    # Generates UMAP and default clustering.
    sc.pp.neighbors(adata, n_pcs=pcs_neighbours)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, key_added=CLUSTER_KEY, resolution=1.0)

    print("\tPlotting data...")
    fig = sc.pl.umap(
        adata,
        color=[
            CLUSTER_KEY,
            REPLICATE_KEY,
            SAMPLE_TYPE_KEY,
            "n_counts",
            "percent_mito",
            "percent_ribo",
            "doublet_score",
        ],
        wspace=1,
        show=False,
        return_fig=True,
    )
    fig.savefig(os.path.join(output_folder_path, "default_umap.svg"))
    plt.close(fig)

    print("\tWriting data to file...")
    adata.write(
        os.path.join(output_folder_path, "clustered.h5ad"),
        compression="gzip",
    )

    print("\tPlotting sample types...")
    replicates = adata.obs[SAMPLE_TYPE_KEY].cat.categories
    for replicate in replicates:
        fig_replicate = sc.pl.umap(
            adata,
            color=[SAMPLE_TYPE_KEY],
            groups=replicate,
            wspace=1,
            show=False,
            return_fig=True,
            title=replicate,
        )
        fig_replicate.savefig(
            os.path.join(
                output_folder_path,
                pathvalidate.sanitize_filename(f"sample_{replicate}.svg"),
            )
        )
        plt.close(fig_replicate)

    print("\tPlotting replicates...")
    replicates = adata.obs[REPLICATE_KEY].cat.categories
    for replicate in replicates:
        fig_replicate = sc.pl.umap(
            adata,
            color=[REPLICATE_KEY],
            groups=replicate,
            wspace=1,
            show=False,
            return_fig=True,
            title=replicate,
        )
        fig_replicate.savefig(
            os.path.join(
                output_folder_path,
                pathvalidate.sanitize_filename(f"replicate_{replicate}.svg"),
            )
        )
        plt.close(fig_replicate)

    print("\tPlotting marker genes...")
    sc.tl.rank_genes_groups(
        adata,
        groupby=CLUSTER_KEY,
        method="wilcoxon",
        key_added="marker_genes_leiden",
    )

    fig = sc.pl.rank_genes_groups_dotplot(
        adata,
        groupby=CLUSTER_KEY,
        standard_scale="var",
        n_genes=5,
        key="marker_genes_leiden",
        show=False,
        return_fig=True,
    )
    fig.savefig(os.path.join(output_folder_path, "default_marker_genes.svg"))
