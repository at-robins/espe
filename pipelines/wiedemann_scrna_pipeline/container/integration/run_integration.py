#!/usr/bin/python
"""This module integrates."""

import anndata
import csv
import json
import logging
import numpy as np
import os
import pathvalidate
import scanpy as sc
import seaborn as sns
import warnings

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER = MOUNT_PATHS["dependencies"]["doublet_detection"] + "/"
SAMPLE_INFO_DIRECTORY = "sample directory"
SAMPLE_INFO_TYPE = "sample type"
DEFAULT_BATCH_KEY = "batch"

# Setup of scanpy.
sc.settings.verbosity = 2
sc.settings.set_figure_params(
    dpi=300,
    facecolor="white",
    frameon=False,
)
sc.settings.figdir = MOUNT_PATHS["output"]


# Sets batch information based on the supplied sample information.
print("Parsing sample information...")
sample_info_path = f"{MOUNT_PATHS['input']}/sample_information.csv"
sample_info_map = {}
with open(sample_info_path, newline="", encoding="utf-8") as csvfile:
    info_reader = csv.DictReader(csvfile, dialect="unix", delimiter=",", quotechar='"')
    for row in info_reader:
        sample_type = sample_info_map.get(row["sample type"])
        sample_directory = row["sample directory"]
        if sample_type is None:
            sample_info_map[row["sample type"]] = [sample_directory]
        else:
            sample_type.append(sample_directory)

for key, values in sample_info_map.items():
    print(f"Processing batch {key} : {values}...", flush=True)
    adatas = []
    for input_directory in values:
        input_path = os.path.join(
            INPUT_FOLDER, input_directory, "filtered_feature_bc_matrix.h5ad"
        )
        print(f"\tReading {input_path}...")
        adatas.append(anndata.read_h5ad(input_path))

    print("\tMerging data...")
    adata_merged = anndata.concat(
        adatas,
        axis=0,
        join="outer",
        merge=None,
        label=DEFAULT_BATCH_KEY,
        keys=values,
    )
    print("\tWriting merged data to file...")
    output_folder_path = os.path.join(
        MOUNT_PATHS["output"], pathvalidate.sanitize_filename(key)
    )
    os.makedirs(output_folder_path, exist_ok=True)
    adata_merged.write(
        f"{output_folder_path}/filtered_feature_bc_matrix.h5ad", compression="gzip"
    )
