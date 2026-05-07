#!/usr/bin/python
"""This module performs custom cell type annotation."""

import anndata

print("Loading decoupler...", flush=True)
import decoupler
import json
import os
import pandas as pd
import scanpy as sc

from matplotlib import pyplot as plt
from matplotlib.colors import TwoSlopeNorm


MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
MARKER_FILE = os.path.join(MOUNT_PATHS["input"], "celltypes.csv")
INPUT_FOLDER = next(iter(MOUNT_PATHS["dependencies"].values()))

# Setup of scanpy.
sc.settings.verbosity = 2
sc.set_figure_params(
    scanpy=True,
    # In case of bitmap exports, use high quality.
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
sc.settings.figdir = MOUNT_PATHS["output"]


def process_data(file_path_input, output_folder_path, network):
    """
    Runs the ULM.
    """
    print(f"Processing file {file_path_input}", flush=True)
    print("\tReading data...", flush=True)
    adata = anndata.read_h5ad(file_path_input)
    print("\tRunning ULM...", flush=True)
    decoupler.mt.ulm(data=adata, net=network, tmin=3)
    score = decoupler.pp.get_obsm(adata, key="score_ulm")

    print("\tPlotting data...", flush=True)
    all_umap_samples = score.var_names
    chunk_size = 4 * 8
    for chunk_index in range(0, len(all_umap_samples), chunk_size):
        fig = sc.pl.umap(
            score,
            color=all_umap_samples[chunk_index : chunk_index + chunk_size],
            wspace=1,
            show=False,
            return_fig=True,
            cmap="RdBu_r",
            vcenter=0,
        )
        fig.savefig(os.path.join(output_folder_path, f"cell_types_{chunk_index}.svg"))
        plt.close(fig)


if not os.path.isfile(MARKER_FILE):
    print("No cell type marker file present. Exiting...", flush=True)
    sys.exit()

print("Loading cell type network...", flush=True)
celltype_network = pd.read_csv(
    MARKER_FILE, sep=",", header=0, index_col=None, encoding="utf-8"
)

# Iterates over all directories and loads sample information.
print("Searching for data...", flush=True)
for root, dirs, files in os.walk(INPUT_FOLDER):
    for file in files:
        if file.casefold().endswith(".h5ad"):
            file_path_input = os.path.join(root, file)
            output_folder_path = os.path.join(
                MOUNT_PATHS["output"],
                os.path.normpath(os.path.relpath(root, INPUT_FOLDER)),
            )
            os.makedirs(output_folder_path, exist_ok=True)
            process_data(file_path_input, output_folder_path, celltype_network)
