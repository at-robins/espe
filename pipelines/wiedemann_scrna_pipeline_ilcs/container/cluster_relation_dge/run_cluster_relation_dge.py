#!/usr/bin/python
"""This module calculates differentially expressed genes."""

import anndata
import anndata2ri
import csv
import functools
import json
import logging
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

from matplotlib import pyplot as plt
from rpy2.robjects.packages import importr

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER_COUNTS = MOUNT_PATHS["dependencies"]["clustering_cell_type"] + "/"
INPUT_FOLDER_SAMPLING = MOUNT_PATHS["dependencies"]["cluster_resolution_data"] + "/"
INPUT_FOLDER_TREE = MOUNT_PATHS["dependencies"]["cluster_relation_tree"] + "/"

# The replicates.
KEY_COUNTS_REPLICATE = "replicate_name"
KEY_PSEUDOBULK_REPLICATE = "replicate"
KEY_PSEUDOBULK_SAMPLE = "sample"
KEY_PSEUDOBULK_N_CELLS = "cellcount"
VALUE_PSEUDOBULK_SAMPLE_TEST = "test"
VALUE_PSEUDOBULK_SAMPLE_REFERENCE = "reference"
VALID_R_CHARACTERS = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_"

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
    vector_friendly=True,
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


def load_sampling(sampling_path, tree):
    """
    Loads the cluster sampling file and appends the info to the tree data.
    """
    print("\tLoading cluster sampling data...", flush=True)
    with open(sampling_path, newline="", encoding="utf-8") as csvfile:
        info_reader = csv.reader(csvfile, dialect="unix", delimiter=",", quotechar='"')
        info_added = 0
        for row in info_reader:
            resolution = float(row[0])
            # Floating point equality works here as we are looking for
            # the exact same resolution and both floats have been parsed
            # from the same string.
            tree_entry = next(
                (entry for entry in tree if entry["resolution"] == resolution), None
            )
            if tree_entry is not None:
                tree_entry["clustering"] = list(map(int, row[1:]))
                info_added += 1

        if info_added != len(tree):
            raise ValueError(
                (
                    f"The cluster sampling file {sampling_path} does not"
                    "contain all resolution entries present in the cluster relation tree."
                )
            )


def load_tree(tree_path):
    """
    Loads the cluster relation tree file.
    """
    print("\tLoading cluster relation tree...", flush=True)
    with open(tree_path, mode="rt", encoding="utf-8") as tree:
        return json.load(tree)


def n_cluster_obs_name(number_of_clusters):
    """
    Returns the observation name for the clustering info based on the number of clusters specified.
    """
    return f"number_of_clusters_{number_of_clusters}"


def load_counts(counts_path, tree, output_folder):
    """
    Loads the count matrix and plots the sampled clusters.
    """
    print("\tLoading count data...", flush=True)
    adata = anndata.read_h5ad(counts_path)
    entry_obs_names = []
    for entry in tree:
        entry_obs_name = n_cluster_obs_name(entry["number_of_clusters"])
        adata.obs[entry_obs_name] = pd.Categorical(entry["clustering"])
        entry_obs_names.append(entry_obs_name)

    print("\tPlotting data...")
    fig = sc.pl.umap(
        adata,
        color=entry_obs_names,
        legend_loc="on data",
        show=False,
        return_fig=True,
    )
    fig.savefig(os.path.join(output_folder, "umap.svg"))
    plt.close(fig)

    adata.X = adata.layers["counts"]
    return adata


def aggregate_pseudobulk_data(adata, mask):
    """
    Aggregates pseudobulk data.
    """
    batches = adata.obs[KEY_COUNTS_REPLICATE].cat.categories
    pseudobulk_dataframe = pd.DataFrame()
    replicate_array = []
    n_obs_array = []

    print(f"\tNumber of total observations: {adata.n_obs}")
    for batch in batches:
        print(f"\tProcessing replicate {batch}")
        batch_mask = adata.obs[KEY_COUNTS_REPLICATE] == batch
        final_mask = np.logical_and(batch_mask, mask)
        adata_subset = adata[final_mask]
        n_obs_subset = adata_subset.n_obs
        print(f"\t\tNumber of subset observations: {n_obs_subset}")
        if n_obs_subset > 0:
            # A single batch / replicate there can only have one sample type.
            aggregated_row = adata_subset.to_df().agg(np.sum)
            replicate_array.append(batch)
            n_obs_array.append(n_obs_subset)
            pseudobulk_dataframe = pd.concat(
                [pseudobulk_dataframe, aggregated_row.to_frame().T],
                ignore_index=False,
                join="outer",
            )
        else:
            print("\t\tNo observations found. Skipping subset...")

    return (
        pseudobulk_dataframe,
        replicate_array,
        n_obs_array,
    )


def dge_for_split(
    adata, sample_cluster_id, reference_cluster_ids, number_of_clusters, output_path
):
    """
    Performs differential gene expression analysis for splits in the cluster tree relations.
    """
    obs_name = n_cluster_obs_name(number_of_clusters)
    test_cluster_mask = adata.obs[obs_name] == sample_cluster_id
    # Combines all the reference cluster IDs
    reference_cluster_mask = np.zeros((len(test_cluster_mask),), dtype=bool)
    for reference_id in reference_cluster_ids:
        reference_cluster_mask = np.logical_or(
            reference_cluster_mask, adata.obs[obs_name] == reference_id
        )
    (
        tmp_pseudobulk_test,
        replicates,
        n_cells,
    ) = aggregate_pseudobulk_data(adata=adata, mask=test_cluster_mask)
    (
        tmp_pseudobulk_reference,
        tmp_replicates_reference,
        tmp_n_cells_reference,
    ) = aggregate_pseudobulk_data(adata=adata, mask=reference_cluster_mask)

    if len(tmp_pseudobulk_test) < 2:
        print(
            f"\tThe test sample has only {len(tmp_pseudobulk_test)} replicates. This is not enough to measure dispersion. Skipping cluster comparison..."
        )
        return

    if len(tmp_pseudobulk_reference) < 2:
        print(
            f"\tThe reference sample has only {len(tmp_pseudobulk_reference)} replicates. This is not enough to measure dispersion. Skipping cluster comparison..."
        )
        return

    pseudobulk_data = pd.concat(
        [tmp_pseudobulk_test, tmp_pseudobulk_reference],
        ignore_index=False,
        join="outer",
    )
    samples = [VALUE_PSEUDOBULK_SAMPLE_TEST] * len(replicates)
    samples.extend([VALUE_PSEUDOBULK_SAMPLE_REFERENCE] * len(tmp_replicates_reference))
    replicates.extend(tmp_replicates_reference)
    n_cells.extend(tmp_n_cells_reference)

    # Replaces NAs with zeros.
    pseudobulk_data.fillna(0, inplace=True)

    # Sets row names.
    row_names = list(
        map(
            convert_string_to_r,
            [f"{sample}_{replicate}" for sample, replicate in zip(samples, replicates)],
        )
    )
    pseudobulk_data.index = row_names

    # Sets additional variables.
    pseudobulk_adata = anndata.AnnData(X=pseudobulk_data)
    pseudobulk_adata.obs[KEY_PSEUDOBULK_SAMPLE] = samples
    pseudobulk_adata.obs[KEY_PSEUDOBULK_REPLICATE] = replicates
    pseudobulk_adata.obs[KEY_PSEUDOBULK_N_CELLS] = n_cells

    # Sets categorical data to a categorical data type.
    pseudobulk_adata.obs[KEY_PSEUDOBULK_SAMPLE] = pseudobulk_adata.obs[
        KEY_PSEUDOBULK_SAMPLE
    ].astype("category")
    pseudobulk_adata.obs[KEY_PSEUDOBULK_REPLICATE] = pseudobulk_adata.obs[
        KEY_PSEUDOBULK_REPLICATE
    ].astype("category")

    output_path_full = os.path.join(
        output_path,
        pathvalidate.sanitize_filename(str(number_of_clusters)),
        pathvalidate.sanitize_filename(
            f"{sample_cluster_id}__vs__{'_'.join(map(str, reference_cluster_ids))}"
        ),
    )
    os.makedirs(output_path_full, exist_ok=True)

    print("\tLoading R dependencies...")
    importr("Seurat")
    importr("edgeR")
    print("\tRunning R...")
    edger_function = ro.r(
        """
        function(data, sample_reference, sample_test, output_path){
           # Creates an edgeR object with counts and grouping factor.
            y <- DGEList(assay(data, "X"), group = colnames(data))
            # Filters out genes with low counts.
            cat("\\tFeatures and samples before subsetting: ", dim(y)[1], " x ", dim(y)[2], "\\n", sep="")
            keep <- filterByExpr(y)
            y <- y[keep, , keep.lib.sizes=FALSE]
            cat("\\tFeatures and samples after subsetting: ", dim(y)[1], " x ", dim(y)[2], "\\n", sep="")
            # Performs normalisation.
            y <- calcNormFactors(y)
            # Creates a vector that is concatentation of condition and variables that we will later use with contrasts
            sample <- colData(data)$sample
            # Creates the design matrix.
            design <- model.matrix(~ 0 + sample)
            # Estimates dispersion.
            y <- estimateDisp(y, design = design)
            # Fits the model.
            fit <- glmQLFit(y, design)
            cat("\\tPlotting data...\\n", sep="")
            svg(paste(output_path, "mds_plot.svg", sep = "/"))
            plotMDS(y)
            dev.off()
            svg(paste(output_path, "bcv_plot.svg", sep = "/"))
            plotBCV(y)
            dev.off()
            # eval workaround as make makeContrasts does not accept a string variable.
            contrast_string = paste("sample", sample_test, " - sample", sample_reference, sep = "")
            cmd <- paste("sample_contrast <- makeContrasts(", contrast_string, ", levels = y$design)", sep ='"')
            eval(parse(text = cmd))
            qlf <- glmQLFTest(fit, contrast=sample_contrast)
            # Returns all of the DE genes and calculates Benjamini-Hochberg adjusted FDR.
            tt <- topTags(qlf, n = Inf)
            tt <- tt$table
            write.csv(tt, paste(output_path, "differential_gene_expression.csv", sep = "/"), row.names=TRUE)
            svg(paste(output_path, "smear_plot.svg", sep = "/"))
            plotSmear(qlf, de.tags = rownames(tt)[which(tt$FDR<=0.05)])
            dev.off()
        }
        """
    )
    edger_function(
        pseudobulk_adata,
        VALUE_PSEUDOBULK_SAMPLE_REFERENCE,
        VALUE_PSEUDOBULK_SAMPLE_TEST,
        output_path_full,
    )


def dge_for_all_splits(adata, tree, output_path):
    """
    Performs differential gene expression analysis for splits in the cluster tree relations.
    """
    for entry_index, entry in enumerate(tree):
        for node in entry["nodes"]:
            child_clusters = node["child_clusters"]
            if len(child_clusters) > 1:
                for child_cluster_index, child_cluster in enumerate(child_clusters):
                    # The reference is all child clusters without the current one.
                    reference = [
                        x
                        for i, x in enumerate(child_clusters)
                        if i != child_cluster_index
                    ]
                    dge_for_split(
                        adata=adata,
                        sample_cluster_id=child_cluster,
                        reference_cluster_ids=reference,
                        number_of_clusters=tree[entry_index + 1]["number_of_clusters"],
                        output_path=output_path,
                    )


print("Searching for data...")
# Iterates over all sample directories and processes them conserving the directory structure.
for root, dirs, files in os.walk(INPUT_FOLDER_COUNTS):
    for file in files:
        if file.casefold().endswith("clustered.h5ad"):
            file_path_counts = os.path.join(root, file)
            file_path_sampling = os.path.join(
                INPUT_FOLDER_SAMPLING,
                root.removeprefix(INPUT_FOLDER_COUNTS),
                "cluster_resolution_data.csv",
            )
            file_path_tree = os.path.join(
                INPUT_FOLDER_TREE,
                root.removeprefix(INPUT_FOLDER_COUNTS),
                "genealogy_cluster_resolution_data.json",
            )
            output_folder_path = os.path.join(
                MOUNT_PATHS["output"], root.removeprefix(INPUT_FOLDER_COUNTS)
            )
            os.makedirs(output_folder_path, exist_ok=True)
            print(
                f"Processing files {file_path_counts}, {file_path_sampling} and {file_path_tree}",
                flush=True,
            )
            relation_tree = load_tree(file_path_tree)
            load_sampling(file_path_sampling, relation_tree)
            counts = load_counts(
                counts_path=file_path_counts,
                tree=relation_tree,
                output_folder=output_folder_path,
            )
            dge_for_all_splits(
                adata=counts, tree=relation_tree, output_path=output_folder_path
            )
