#!/usr/bin/python
"""This module predicts cell-cell communication."""

import anndata
import anndata2ri
import csv
import json
import logging
import math
import multiprocessing
import numpy as np
import os
import pandas as pd
import pathvalidate
import rpy2.rinterface_lib.callbacks as rcb
import rpy2.robjects as ro
import scanpy as sc
import seaborn as sns
import shutil
import warnings

from rpy2.robjects.packages import importr

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER = next(iter(MOUNT_PATHS["dependencies"].values()))
VALID_R_CHARACTERS = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_"

# Loads environment /system variables.
ENV_ORGANISM = os.environ.get("GLOBAL_ORGANISM")
THREADS = math.floor(multiprocessing.cpu_count() * 0.8)
if THREADS < 1:
    THREADS = 1

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


def convert_string_to_r(val: str) -> str:
    """
    Converts a string to a valid R string.
    """
    return "".join(
        list(
            map(
                lambda letter: letter if letter in VALID_R_CHARACTERS else ".",
                val.replace(" ", "_"),
            )
        )
    )


def is_cluster_obs(obs_name):
    """
    Returns True if the observation name contains cluster information.
    """
    return obs_name.startswith("number_of_clusters_")


def cluster_name_to_cluster(obs_name):
    """
    Returns the cluster number from the observation name.
    """
    return obs_name.removeprefix("number_of_clusters_")


def prepare_cluster_name_for_R(cluster_number):
    """
    Returns a valid R / CellChat cluster name.
    This is needed because CellChat does not support '0' as cluster name.
    """
    return f"C{cluster_number}"


def run_cell_communication(
    adata,
    output_path,
):
    """
    Performes the cell-cell communication prediction.
    """

    print("\tLoading R dependencies...")
    importr("CellChat")
    importr("patchwork")
    print("\tRunning R...")
    cell_chat_function = ro.r(
        """
        function(data_matrix, cells, genes, clusters, organism, output_path, n_threads){
            # CellChat prerequisites.
            options(stringsAsFactors = FALSE)
            set.seed(42)
            rownames(data_matrix) <- genes
            colnames(data_matrix) <- cells
            if (future::supportsMulticore()) {
                cat("\\tUsing ", n_threads, " worker threads.", "\\n", sep="")
                future::plan(future::multicore, workers = n_threads)
            } else {
                future::plan(future::multisession)
            }

            # Creates the CellChat object.
            meta = data.frame(labels = clusters, row.names = colnames(data_matrix))
            cellchat <- createCellChat(object = data_matrix, meta = meta, group.by = "labels")

            # Loads the database.
            if (organism == "human") {
                cat("\\tLoading human database...", "\\n", sep="")
                CellChatDB <- CellChatDB.human
            } else {
                cat("\\tLoading mouse database...", "\\n", sep="")
                CellChatDB <- CellChatDB.mouse
            }
            # Subsets to only use secreted signalling.
            cat("\\tSubsetting database...", "\\n", sep="")
            cellchat@DB <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation")

            # Subsets relevant genes.
            cat("\\tFiltering relevant genes...", "\\n", sep="")
            cellchat <- subsetData(cellchat)

            cat("\\tIdentifying overexpressed genes...", "\\n", sep="")
            cellchat <- identifyOverExpressedGenes(cellchat)

            cat("\\tIdentifying overexpressed interactions...", "\\n", sep="")
            cellchat <- identifyOverExpressedInteractions(cellchat)

            cat("\\tComputing communication probability...", "\\n", sep="")
            cellchat <- computeCommunProb(cellchat, type = "triMean")

            cat("\\tFiltering out low abundance communities...", "\\n", sep="")
            cellchat <- filterCommunication(cellchat, min.cells = 10)

            cat("\\tSaving communication probability data...", "\\n", sep="")
            write.csv(
                subsetCommunication(cellchat),
                paste(output_path, "communication_probability.csv", sep = "/"),
                row.names=TRUE
            )

            cat("\\tComputing pathway communication probability...", "\\n", sep="")
            cellchat <- computeCommunProbPathway(cellchat)

            cat("\\tSaving communication probability data...", "\\n", sep="")
            write.csv(subsetCommunication(
                cellchat, slot.name = "netP"),
                paste(output_path, "communication_probability_pathways.csv", sep = "/"),
                row.names=TRUE
            )

            cat("\\tAggregating communication network...", "\\n", sep="")
            cellchat <- aggregateNet(cellchat)

            groupSize <- as.numeric(table(cellchat@idents))
            svg(paste(output_path, "aggregated_network.svg", sep = "/"))
            par(mfrow = c(1,2), xpd=TRUE)
            netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
            netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
            dev.off()
        }
        """
    )

    cells = list(map(str, adata.obs_names))
    genes = list(map(str, adata.var_names))
    for obs_name in adata.obs:
        if is_cluster_obs(obs_name=obs_name):
            cluster = cluster_name_to_cluster(obs_name=obs_name)
            # With a single cluster there cannot be any communication.
            if int(cluster) > 1:
                print(f"\tProcessing {cluster} clusters...", flush=True)
                cluster_labels = list(
                    map(prepare_cluster_name_for_R, adata.obs[obs_name])
                )
                final_output_folder = os.path.join(output_path, cluster)
                os.makedirs(final_output_folder, exist_ok=True)
                cell_chat_function(
                    adata.layers["log1p_norm"].T,
                    ro.StrVector(cells),
                    ro.StrVector(genes),
                    ro.StrVector(cluster_labels),
                    ENV_ORGANISM,
                    final_output_folder,
                    THREADS,
                )


print("Searching for data...")
# Iterates over all sample directories and processes them conserving the directory structure.
for root, dirs, files in os.walk(INPUT_FOLDER):
    for file in files:
        if file.casefold().endswith("cluster_relation.h5ad"):
            file_path_counts = os.path.join(root, file)
            output_folder_path = os.path.join(
                MOUNT_PATHS["output"],
                os.path.normpath(os.path.relpath(root, INPUT_FOLDER)),
            )
            os.makedirs(output_folder_path, exist_ok=True)
            print(
                f"Processing file {file_path_counts}...",
                flush=True,
            )
            adata_clustered = anndata.read_h5ad(file_path_counts)
            run_cell_communication(
                adata=adata_clustered, output_path=output_folder_path
            )
