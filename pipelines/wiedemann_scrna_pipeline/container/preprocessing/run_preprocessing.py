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

sc.settings.verbosity = 2
sc.settings.set_figure_params(
    dpi=80,
    facecolor="white",
    frameon=False,
)

adata = sc.read_10x_h5(filename=f"{INPUT_FOLDER}/sample_filtered_feature_bc_matrix.h5")
adata

adata.var_names_make_unique()
adata

# mitochondrial genes
adata.var["mt"] = adata.var_names.str.lower().startswith("mt-")
# ribosomal genes
adata.var["ribo"] = adata.var_names.str.lower().startswith(("rps", "rpl"))
# hemoglobin genes.
adata.var["hb"] = adata.var_names.str.contains(("^HB[^(P)]"))

sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True
)
adata

p1 = sns.displot(adata.obs["total_counts"], bins=100, kde=False)
# sc.pl.violin(adata, 'total_counts')
p2 = sc.pl.violin(adata, "pct_counts_mt")
p3 = sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")

p1.savefig(f"{MOUNT_PATHS['output']}/p1.png")
p2.savefig(f"{MOUNT_PATHS['output']}/p2.png")
p3.savefig(f"{MOUNT_PATHS['output']}/p3.png")

