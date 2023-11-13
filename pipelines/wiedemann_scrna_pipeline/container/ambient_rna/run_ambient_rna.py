#!/usr/bin/python
"""This module removes ambient RNA."""

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

from rpy2.robjects.packages import importr

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER = MOUNT_PATHS["dependencies"]["preprocessing"] + "/"

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
    dpi=80,
    facecolor="white",
    frameon=False,
)
sc.settings.figdir = MOUNT_PATHS["output"]


def get_raw_file(raw_file_folder):
    """
    Returns the raw data file in the specified folder or None if
    none is present.

    Keyword arguments:
    raw_file_folder -- the directory path to search the raw file in
    """
    for file_path in os.listdir(raw_file_folder):
        full_path = os.path.join(raw_file_folder, file_path)
        if os.path.isfile(full_path) and file_path.casefold().endswith(
            "raw_feature_bc_matrix.h5"
        ):
            return full_path
    return None


def process_data(file_path_filtered, file_path_raw, output_folder_path, skip_soupx, metrics_writer):
    """
    Removes ambient RNA.
    """
    print(f"Processing files {file_path_filtered} and {file_path_raw}", flush=True)
    print("\tReading filtered data...")
    adata_filtered = anndata.read_h5ad(file_path_filtered)
    if skip_soupx:
        print("\tSkipping SoupX...")
    else:
        gene_ids_filtered = adata_filtered.var["gene_ids"].to_numpy()

        print("\tNormalising data...")
        adata_filtered_tmp = adata_filtered.copy()
        sc.pp.normalize_per_cell(adata_filtered_tmp)
        sc.pp.log1p(adata_filtered_tmp)

        print("\tClustering...")
        sc.pp.pca(adata_filtered_tmp)
        sc.pp.neighbors(adata_filtered_tmp)
        sc.tl.leiden(adata_filtered_tmp, key_added="soupx_groups")
        soupx_groups = adata_filtered_tmp.obs["soupx_groups"]
        # Deletes the data reference to save memory.
        del adata_filtered_tmp
        soupx_cells = adata_filtered.obs_names
        soupx_genes = adata_filtered.var_names
        soupx_data_filtered = adata_filtered.X.T

        print("\tReading and preprocessing raw data...")
        adata_raw = sc.read_10x_h5(file_path_raw)
        gene_ids_raw = adata_raw.var["gene_ids"].to_numpy()
        if gene_ids_filtered.size < gene_ids_raw.size:
            warnings.warn(
                (
                    f"The filtered feature matrix ({gene_ids_filtered.size}) "
                    f"has less features than the raw matrix ({gene_ids_raw.size}). "
                    "Trying to filter the raw matrix..."
                )
            )
            raw_mask = np.in1d(gene_ids_raw, gene_ids_filtered)
            adata_raw = adata_raw[:, raw_mask].copy()
            gene_ids_raw = adata_raw.var["gene_ids"].to_numpy()

        if gene_ids_filtered.size != gene_ids_raw.size:
            raise ValueError(
                (
                    f"The number of features of the filtered matrix ({gene_ids_filtered.size}) "
                    f"does not match the raw matrix ({gene_ids_raw.size}) feature number."
                )
            )

        adata_raw.var_names_make_unique()
        soupx_data_raw = adata_raw.X.T
        # Deletes the data reference to save memory.
        del adata_raw

        print("\tLoading SoupX...")
        importr("SoupX")
        print("\tRunning SoupX...")
        soupx_function = ro.r(
            """
            function(data, data_tod, genes, cells, soupx_groups, output_path) {
                # Constructs dataframes and converts them to sparse matrices.
                rownames(data) = genes
                colnames(data) = cells
                data <- as(data, "sparseMatrix")
                data_tod <- as(data_tod, "sparseMatrix")

                # Generates SoupChannel and sets an additional metadata profile as well as clustering information.
                sc = SoupChannel(data_tod, data, calcSoupProfile = FALSE)
                soupProf = data.frame(
                    row.names = rownames(data),
                    est = rowSums(data)/sum(data),
                    counts = rowSums(data)
                )
                sc = setSoupProfile(sc, soupProf)
                sc = setClusters(sc, soupx_groups)

                # Estimates and plots contamination fraction.
                svg(paste(output_path, "contamination_fraction_plot.svg", sep = "/"))
                sc  = autoEstCont(sc, doPlot=TRUE)
                dev.off()
                # Return corrected, rounded counts.
                return(adjustCounts(sc, roundToInt = TRUE))
            }
            """
        )
        soupx_output = soupx_function(
            soupx_data_filtered,
            soupx_data_raw,
            soupx_genes,
            soupx_cells,
            soupx_groups,
            output_folder_path,
        )
        adata_filtered.layers["counts"] = adata_filtered.X
        adata_filtered.layers["soupX_counts"] = soupx_output.T
        adata_filtered.X = adata_filtered.layers["soupX_counts"]
    n_cells_before_filter = adata_filtered.n_vars
    print(f"Total number of features before filtering: {n_cells_before_filter}")
    sc.pp.filter_genes(adata_filtered, min_cells=20)
    n_cells_after_filter = adata_filtered.n_vars
    print(f"Number of features after filtering: {n_cells_after_filter}")
    print("\tWriting metrics to file...")
    metrics_writer.writerow(
        [file_path_filtered, n_cells_before_filter, n_cells_after_filter]
    )
    print("\tWriting filtered data to file...")
    adata_filtered.write(f"{output_folder_path}/corrected.h5ad", compression="gzip")

with open(
    f"{MOUNT_PATHS['output']}/metrics.csv", mode="w", newline="", encoding="utf-8"
) as csvfile:
    metrics_writer = csv.writer(
        csvfile, dialect="unix", delimiter=",", quotechar='"', quoting=csv.QUOTE_MINIMAL
    )
    metrics_writer.writerow(
        [
            "Sample",
            "Number of features before filtering",
            "Number of features after filtering",
        ]
    )
    # Iterates over all sample directories and processes them conserving the directory structure.
    skip_env = os.environ.get("SKIP")
    skip = skip_env is not None and skip_env == "true"

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
                process_data(
                    file_path_filtered,
                    get_raw_file(folder_path_raw),
                    output_folder_path,
                    skip,
                    metrics_writer
                )
