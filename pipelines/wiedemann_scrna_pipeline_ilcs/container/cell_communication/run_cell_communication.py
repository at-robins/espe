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

# Loads environment /system variables.
ENV_ORGANISM = os.environ.get("GLOBAL_ORGANISM")
THREADS = math.floor(multiprocessing.cpu_count() * 0.8)
if THREADS < 1:
    THREADS = 1
elif THREADS > 4:
    # Limits the number of concurrent threads as the
    # R code that is executed in parallel runs relatively
    # fast even on a low number of threads, while execssive
    # parallelisation induces a hugh memory overhead
    # that might get the process killed.
    THREADS = 4

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
    return f"C{cluster_number:0>2}"


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
    importr("ggplot2")
    importr("ComplexHeatmap")
    importr("NMF")
    importr("ggalluvial")
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
            cluster_count = length(unique(clusters))
            if (cluster_count < 2) {
                cat(
                    "\\tNot enough cell populations (",
                    cluster_count,
                    ") to calculate communication. Skipping the sample...",
                    "\\n",
                    sep=""
                )
                return()
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
            cellchat@DB <- subsetDB(CellChatDB)

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

            cat("\\tComputing centrality score...", "\\n", sep="")
            cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
            
            cat("\\tPlotting...", "\\n", sep="")

            # Accesses all the signaling pathways showing significant communications.
            all_pathways <- cellchat@netP$pathways
            for (pathway in all_pathways) {
                # Plots the interaction probability heatmap.
                svg(
                    paste(output_path, "/heatmap_probability_", pathway, ".svg", sep = ""),
                    width = 0.28 * cluster_count,
                    height = 0.28 * cluster_count,
                )
                draw(netVisual_heatmap(cellchat, signaling = pathway, color.heatmap = "Reds"))
                dev.off()

                # Plots the chord diagram.
                svg(paste(output_path, "/chord_", pathway, ".svg", sep = ""))
                netVisual_aggregate(cellchat, signaling = pathway, layout = "chord")
                dev.off()

                # Computes and visualize the contribution of 
                # each ligand-receptor pair to the overall signaling pathway.
                contribution_plot <- netAnalysis_contribution(cellchat, signaling = pathway)
                ggsave(
                    filename=paste(output_path, "/contribution_", pathway, ".svg", sep = ""),
                    plot=contribution_plot
                )

                # Plots centrality score.
                svg(
                    paste(output_path, "/centrality_", pathway, ".svg", sep = ""),
                    width = 1.1 + 0.15 * cluster_count,
                    height = 1.8
                )
                # Function uses metric system internally.
                netAnalysis_signalingRole_network(
                    cellchat,
                    signaling = pathway,
                    width = 0.15 * cluster_count * 2.54,
                    height = 1.0 * 2.54
                )
                dev.off()
            }

            # Plots signal contribution
            svg(
                paste(output_path, "/heatmap_contribution_outgoing.svg", sep = ""),
                width = (1.6 + 0.235 * cluster_count) * 2.54
            )
            # Function uses metric system internally.
            draw(netAnalysis_signalingRole_heatmap(
                cellchat,
                pattern = "outgoing",
                width = (1.6 + 0.235 * cluster_count) * 2.54
            ))
            dev.off()
            svg(
                paste(output_path, "/heatmap_contribution_incoming.svg", sep = ""),
                width = (1.6 + 0.235 * cluster_count) * 2.54
            )
            # Function uses metric system internally.
            draw(netAnalysis_signalingRole_heatmap(
                cellchat,
                pattern = "incoming",
                width = (1.6 + 0.235 * cluster_count) * 2.54
            ))
            dev.off()

            cat("\\tInferring patterns...", "\\n", sep="")
            pattern_min = 2
            pattern_max = 10
            # Reduces the number of inferred patterns if uninformative.
            if (pattern_max >= cluster_count) {
                if (cluster_count > pattern_min) {
                    pattern_max = cluster_count - 1
                } else {
                    pattern_max = pattern_min
                }
            }
            outgoing_pattern_plot = selectK(
                cellchat,
                k.range = seq(pattern_min, pattern_max),
                pattern = "outgoing"
            )
            ggsave(
                filename=paste(output_path, "/outgoing_pattern_inference.svg", sep = ""),
                plot=outgoing_pattern_plot
            )  
            incoming_pattern_plot = selectK(
                cellchat,
                k.range = seq(pattern_min, pattern_max),
                pattern = "incoming"
            )
            ggsave(
                filename=paste(output_path, "/incoming_pattern_inference.svg", sep = ""),
                plot=incoming_pattern_plot
            )
            for (pattern_count in pattern_min:pattern_max) {
                cat("\\tPlotting ", pattern_count, " patterns...", "\\n", sep="")
                tryCatch({
                    svg(
                        paste(output_path, "/communication_pattern_heatmap_outgoing_", pattern_count, ".svg", sep = ""),
                        height = 1.3 + 0.3 * cluster_count
                    )
                    # The function internally plots two plots next to
                    # each other using metric units.
                    cellchat <- identifyCommunicationPatterns(
                        cellchat,
                        pattern = "outgoing",
                        k = pattern_count,
                        height = 0.14 * cluster_count * 2.54
                    )
                    dev.off()
                    ggsave(
                        filename=paste(
                            output_path,
                            "/communication_pattern_river_outgoing_",
                            pattern_count,
                            ".svg",
                            sep = ""
                        ),
                        plot=netAnalysis_river(cellchat, pattern = "outgoing"),
                        width = 7,
                        height = 1.4 + 0.224 * cluster_count
                    )
                    ggsave(
                        filename=paste(
                            output_path,
                            "/communication_pattern_dot_outgoing_",
                            pattern_count,
                            ".svg",
                            sep = ""
                        ),
                        plot=netAnalysis_dot(cellchat, pattern = "outgoing")
                    )

                    svg(
                        paste(output_path, "/communication_pattern_heatmap_incoming_", pattern_count, ".svg", sep = ""),
                        height = 1.3 + 0.3 * cluster_count
                    )
                    # The function internally plots two plots next to
                    # each other using metric units.
                    cellchat <- identifyCommunicationPatterns(
                        cellchat,
                        pattern = "incoming",
                        k = pattern_count,
                        height = 0.14 * cluster_count * 2.54
                    )
                    dev.off()
                    ggsave(
                        filename=paste(
                            output_path,
                            "/communication_pattern_river_incoming_",
                            pattern_count,
                            ".svg",
                            sep = ""
                        ),
                        plot=netAnalysis_river(cellchat, pattern = "incoming"),
                        width = 7,
                        height = 1.4 + 0.224 * cluster_count
                    )
                    ggsave(
                        filename=paste(
                            output_path,
                            "/communication_pattern_dot_incoming_",
                            pattern_count,
                            ".svg",
                            sep = ""
                        ),
                        plot=netAnalysis_dot(cellchat, pattern = "incoming")
                    )
                }, error = function(e) {
                    # The pattern count might be too high for the amount of clusters present,
                    # thus execution is skipped and the error logged.
                    cat("\\t\\tThe pattern count is too high. Ignoring expected error.", "\\n", sep="")
                    message(e, "\\n")
                    cat("\\t\\tClosing open graphics devices after error...", "\\n", sep="")
                    while (dev.cur() > 1) {
                        dev.off()
                    }
                })
            }
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
