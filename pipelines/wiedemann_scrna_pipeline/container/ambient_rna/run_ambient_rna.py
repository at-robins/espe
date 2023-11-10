#!/usr/bin/python
"""This module removes ambient RNA."""

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

rcb.logger.setLevel(logging.INFO)
ro.pandas2ri.activate()
anndata2ri.activate()

%load_ext rpy2.ipython

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER = MOUNT_PATHS["dependencies"]["preprocessing"] + "/"

# Setup of scanpy.
sc.settings.verbosity = 2
sc.settings.set_figure_params(
    dpi=80,
    facecolor="white",
    frameon=False,
)
sc.settings.figdir = MOUNT_PATHS["output"]

def get_raw_file(raw_file_folder):
    # TODO: load raw file


def process_data(file_path_filtered, file_path_raw, output_folder_path, metrics_writer):
    """
    Removes ambient RNA.
    """
    print(f"Processing files {file_path_filtered} and {file_path_raw}")
    print("\tReading data...")
    adata = sc.read_10x_h5(filename=file_path)
    print("\tMaking variable names unique...")
    adata.var_names_make_unique()


with open(f"{MOUNT_PATHS['output']}/metrics.csv", mode="w", newline="") as csvfile:
    metrics_writer = csv.writer(
        csvfile, dialect="unix", delimiter=",", quotechar='"', quoting=csv.QUOTE_MINIMAL
    )
    metrics_writer.writerow(
        [
            "Sample",
            "Number of cells before filtering",
            "Number of cells after filtering",
            "Number of QC outlier cells",
            "Number of mitochondiral count outlier cells",
        ]
    )
    # Iterates over all sample directories and processes them conserving the directory structure.
    for root, dirs, files in os.walk(INPUT_FOLDER):
        for file in files:
            if file.casefold().endswith("preprocessed.h5ad"):
                file_path_filtered = os.path.join(root, file)
                folder_path_raw = os.path.join(
                    MOUNT_PATHS["input"], root.removeprefix(INPUT_FOLDER)
                )
                output_folder_path = os.path.join(
                    MOUNT_PATHS["output"], root.removeprefix(INPUT_FOLDER)
                )
                os.makedirs(output_folder_path, exist_ok=True)
                process_data(file_path_filtered, output_folder_path, metrics_writer)
