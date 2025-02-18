#!/usr/bin/python
"""This module runs data preprocessing."""

import csv
import ddqc
import json
import matplotlib
import matplotlib.pyplot as plt
import os
import pegasus as pg
import pegasusio as pgio
import scanpy as sc
from scipy import sparse


MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER = MOUNT_PATHS["input"] + "/"
MITOCHONDRIAL_BASE_PREFIX = "MT-"

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
sc.settings.figdir = ""


def load_anndata(file_path):
    """
    Loads the data.
    """
    if file_path.casefold().endswith("h5"):
        return sc.read_10x_h5(filename=file_path)
    else:
        return sc.read_h5ad(filename=file_path)


def process_data(file_path, output_folder_path, metrics_writer):
    """
    Processes datasets and removes low qulaity reads.

    Keyword arguments:
    file_path -- the input path of the dataset hd5 file
    output_folder_path -- the corresponding output path
    metrics_writer -- the CSV writer to persist calculated metrics
    """
    print(f"Processing file {file_path}", flush=True)
    print("\tReading data...")
    adata = load_anndata(file_path=file_path)
    print("\tMaking variable names unique...")
    adata.var_names_make_unique()
    print("\tEnsuring the matrix is sparse...")
    adata.X = sparse.csr_matrix(adata.X)
    n_cells_raw = adata.n_obs

    print("\tConverting to MultimodalData...", flush=True)
    mmdata = pgio.MultimodalData(adata)
    print("\tDetermining mitochondrial prefix...", flush=True)
    # The casing of the mitochondrial gene prefix is undefined
    # and the QC function does not allow regex, so we use the first
    # gene hit to intrapolate the prefix or use the default if none is found.
    mitochondiral_prefix = MITOCHONDRIAL_BASE_PREFIX
    for gene_name in adata.var_names:
        if gene_name.casefold().startswith(mitochondiral_prefix.casefold()):
            mitochondiral_prefix = gene_name[0 : len(MITOCHONDRIAL_BASE_PREFIX)]
            print(f'\tDetected mitochondrial prefix: "{mitochondiral_prefix}"')
            break
    print("\tRunning QC...", flush=True)
    qc_df = ddqc.ddqc_metrics(
        mmdata,
        clustering_method="leiden",
        method="mad",
        return_df_qc=True,
        random_state=42,
        mito_prefix=mitochondiral_prefix,
    )
    print("\tExporting QC plots...", flush=True)
    plt.savefig(os.path.join(output_folder_path, "ddqc_line_plot_relative.svg"))
    plt.close()
    plt.savefig(os.path.join(output_folder_path, "ddqc_line_plot_absolute.svg"))
    plt.close()
    plt.savefig(
        os.path.join(output_folder_path, "ddqc_box_plot_mitochondrial_reads.svg")
    )
    plt.close()
    plt.savefig(os.path.join(output_folder_path, "ddqc_box_plot_genes.svg"))
    plt.close()
    print("\tExporting QC table...", flush=True)
    qc_df.to_csv(
        os.path.join(output_folder_path, "ddqc_table.csv"),
        sep=",",
        encoding="utf-8",
    )
    print("\tFiltering QC outliers...")
    print(f"\tNumber of cells before filtering: {n_cells_raw}", flush=True)
    pg.filter_data(mmdata)
    print("\tConverting back to AnnData...", flush=True)
    adata = mmdata.to_anndata()
    n_cells_filtered = adata.n_obs
    print(f"\tNumber of cells after filtering: {n_cells_filtered}")

    print("\tWriting metrics to file...", flush=True)
    metrics_writer.writerow(
        [
            file_path.removeprefix(INPUT_FOLDER),
            n_cells_raw,
            n_cells_filtered,
            len(qc_df[qc_df["percent_mito_passed_qc"] == False].index),
            len(qc_df[qc_df["percent_ribo_passed_qc"] == False].index),
            len(qc_df[qc_df["n_counts_passed_qc"] == False].index),
            len(qc_df[qc_df["n_genes_passed_qc"] == False].index),
        ]
    )

    print("\tWriting filtered data to file...")
    # Backup the raw counts in a separate layer.
    adata.layers["counts"] = adata.X
    adata.write(
        os.path.join(output_folder_path, "filtered_feature_bc_matrix.h5ad"),
        compression="gzip",
    )


with open(
    os.path.join(MOUNT_PATHS["output"], "metrics.csv"),
    mode="w",
    newline="",
    encoding="utf-8",
) as csvfile:
    metrics_writer = csv.writer(
        csvfile, dialect="unix", delimiter=",", quotechar='"', quoting=csv.QUOTE_MINIMAL
    )
    metrics_writer.writerow(
        [
            "Sample",
            "Number of cells before filtering",
            "Number of cells after filtering",
            "Number of mitochondiral gene count outlier cells",
            "Number of ribosomal gene count outlier cells",
            "Number of total count outlier cells",
            "Number of gene variance outlier cells",
        ]
    )
    # Iterates over all sample directories and processes them conserving the directory structure.
    for root, dirs, files in os.walk(INPUT_FOLDER):
        for file in files:
            if "filtered_feature_bc_matrix.h5" in file.casefold():
                input_file_path = os.path.join(root, file)
                output_folder_path = os.path.join(
                    MOUNT_PATHS["output"], root.removeprefix(INPUT_FOLDER)
                )
                os.makedirs(output_folder_path, exist_ok=True)
                process_data(input_file_path, output_folder_path, metrics_writer)
