#!/usr/bin/python
"""This module integrates."""

import anndata
import csv
import json
import logging
import numpy as np
import os
import pandas as pd
import pathvalidate
import scanpy as sc
import seaborn as sns
import warnings

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER = next(iter(MOUNT_PATHS["dependencies"].values()))
SAMPLE_INFO_DIRECTORY = "sample directory"
SAMPLE_INFO_TYPE = "sample type"
REPLICATE_KEY = "replicate_name"
SAMPLE_TYPE_KEY = "sample_type"

# Setup of scanpy.
sc.settings.verbosity = 2

# Sets batch information based on the supplied sample information.
print("Parsing sample information...")
sample_info_path = f"{MOUNT_PATHS['input']}/sample_information.csv"
adatas = []
with open(sample_info_path, newline="", encoding="utf-8") as csvfile:
    info_reader = csv.DictReader(csvfile, dialect="unix", delimiter=",", quotechar='"')
    for row in info_reader:
        sample_type = row["sample type"]
        sample_directory = row["sample directory"]
        input_path = os.path.join(
            INPUT_FOLDER, sample_directory, "filtered_feature_bc_matrix.h5ad"
        )
        print(f"\tReading {input_path}...")
        sample_data = anndata.read_h5ad(input_path)
        sample_data.obs[REPLICATE_KEY] = pd.Categorical(
            np.repeat(sample_directory, sample_data.n_obs)
        )
        sample_data.obs[SAMPLE_TYPE_KEY] = pd.Categorical(
            np.repeat(sample_type, sample_data.n_obs)
        )
        adatas.append(sample_data)

print("\tMerging data...")
adata_merged = anndata.concat(
    adatas,
    axis=0,
    join="outer",
    merge="same",
    uns_merge="same",
    label=None,
    keys=None,
)
adata_merged.obs_names_make_unique()

print("\tWriting merged data to file...")
output_file_path = os.path.join(MOUNT_PATHS["output"], "merged.h5ad")
os.makedirs(MOUNT_PATHS["output"], exist_ok=True)
adata_merged.write(output_file_path, compression="gzip")
