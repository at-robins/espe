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
INPUT_FOLDER = next(iter(MOUNT_PATHS["dependencies"].values()))
TOP_FEATURES = 4000
REPLICATE_NAME_KEY = "replicate_name"
SAMPLE_TYPE_KEY = "sample_type"

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
        REPLICATE_NAME_KEY in adata.obs_keys()
        and adata.obs[REPLICATE_NAME_KEY].cat.categories.size > 1
    )


def calculate_deviance_unbatched(adata: anndata.AnnData):
    """
    Returns binomial deviance.
    """
    print("\t\tLoading R dependencies...")
    importr("scry")
    print("\t\tRunning R code...")
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
    return np.nan_to_num(scry_function(adata).T)


def get_deviance_mask_unbatched(adata: anndata.AnnData):
    """
    Returns the mask for the top deviant genes.
    """
    binomial_deviance = calculate_deviance_unbatched(adata)
    sorted_indices = binomial_deviance.argsort()
    if sorted_indices.size > TOP_FEATURES:
        print(f"\t\t{sorted_indices.size} features present. Using top {TOP_FEATURES}.", flush=True)
        sorted_indices = sorted_indices[-TOP_FEATURES:]
    else:
        print(f"\t\tOnly {sorted_indices.size} features present. Using all of them.", flush=True)
    mask = np.zeros(adata.var_names.shape, dtype=bool)
    mask[sorted_indices] = True
    return mask


def get_deviance_mask_batched(adata: anndata.AnnData):
    """
    Returns the mask for the top deviant genes for batched data.
    """
    sample_types = adata.obs[SAMPLE_TYPE_KEY].cat.categories
    mask_total = np.zeros(adata.var_names.shape, dtype=bool)
    for sample_type in sample_types:
        print(f"\tCalculating deviance for batch {sample_type}...", flush=True)
        mask_sample_type = None
        sample_type_subset = adata[adata.obs[SAMPLE_TYPE_KEY] == sample_type]
        replicates = sample_type_subset.obs[REPLICATE_NAME_KEY].cat.categories
        for replicate in replicates:
            print(f"\t\tCalculating deviance for replicate {replicate}...", flush=True)
            replicate_subset = sample_type_subset[sample_type_subset.obs[SAMPLE_TYPE_KEY] == sample_type]
            if mask_sample_type is None:
                mask_sample_type = get_deviance_mask_unbatched(replicate_subset)
            else:
                mask_sample_type = np.logical_and(mask_sample_type, get_deviance_mask_unbatched(replicate_subset))
        if mask_sample_type is not None:
            mask_total = np.logical_or(mask_total, mask_sample_type)
    return mask_total


file_path_input = os.path.join(INPUT_FOLDER, "normalised.h5ad")
print(f"Processing file {file_path_input}")
print("\tReading filtered data...", flush=True)
adata = anndata.read_h5ad(file_path_input)
adata.X = adata.layers["counts"]

print("\tCalculating highly deviant genes...")
if is_batch(adata):
    print("\tFound batch information.", flush=True)
    batch_key = REPLICATE_NAME_KEY
    adata.var["highly_deviant"] = get_deviance_mask_batched(adata)
else:
    print("\tNo batch information found. Proceeding without...", flush=True)
    batch_key = None
    adata.var["highly_deviant"] = get_deviance_mask_unbatched(adata)

print("\tCalculating highly variable genes...", flush=True)
sc.pp.highly_variable_genes(adata, layer="scran_normalisation", batch_key=batch_key)

print("\tPlotting data...", flush=True)
fig, ax = plt.subplots()
sns.scatterplot(
    data=adata.var, x="means", y="dispersions", hue="highly_deviant", s=5, ax=ax
)
ax.set_xlim(None, 4)
ax.set_ylim(None, 4)
ax.set(xlabel="Means", ylabel="Dispersion")
legend = ax.get_legend()
# Check prevents index out of bound errors if any one category is absent.
legend_label_count = len(legend.get_texts())
if legend_label_count == 2:
    legend.get_texts()[0].set_text("No")
    legend.get_texts()[1].set_text("Yes")
legend.set_title("Highly deviant")
fig.tight_layout()
os.makedirs(MOUNT_PATHS["output"], exist_ok=True)
fig.savefig(os.path.join(MOUNT_PATHS["output"], "scatter_dispersion.svg"))

print("\tWriting filtered data to file...", flush=True)
adata.write(
    os.path.join(MOUNT_PATHS["output"], "feature_selection.h5ad"),
    compression="gzip",
)
