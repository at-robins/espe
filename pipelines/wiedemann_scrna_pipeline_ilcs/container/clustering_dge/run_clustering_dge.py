#!/usr/bin/python
"""This module performs data clustering."""

import anndata
import csv
import json
import logging
import numpy as np
import os
import pathvalidate
import scanpy as sc
import seaborn as sns

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER_CLUSTERING = MOUNT_PATHS["dependencies"]["clustering_cell_type"] + "/"
INPUT_FOLDER_DGE = MOUNT_PATHS["dependencies"]["differential_gene_expression"] + "/"
CELL_TYPE_KEY = "cell_type"
LOG2_FOLD_CHANGE_CUTOFF = 0.5
FDR_CUTOFF = 0.05

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


def plot_genes(input_directory, output_folder_path, sample, cell_type):
    """
    Plots the genes.
    """
    adata_path = os.path.join(
        INPUT_FOLDER_CLUSTERING,
        pathvalidate.sanitize_filename(sample),
        f"clustered_{pathvalidate.sanitize_filename(cell_type)}.h5ad",
    )
    print(f"\tReading clustering data {adata_path}...")
    if os.path.exists(adata_path):
        adata = anndata.read_h5ad(adata_path)

        print("\tReading differential gene expression data...")
        genes_up = []
        genes_down = []
        dge_path = os.path.join(input_directory, "differential_gene_expression.csv")
        with open(dge_path, newline="", encoding="utf-8") as csvfile:
            dge_reader = csv.DictReader(
                csvfile, dialect="unix", delimiter=",", quotechar='"'
            )
            for row in  dge_reader:
                lfc = float(row["logFC"])
                fdr = float(row["FDR"])
                if np.absolute(lfc) >= LOG2_FOLD_CHANGE_CUTOFF and fdr <= FDR_CUTOFF:
                    gene = row[""] # Genes do not have a header.
                    if lfc > 0:
                        genes_up.append(gene)
                    else:
                        genes_down.append(gene)
   
                    if gene not in adata.var_names:
                        print(f"\tGene {gene} not present in data. Adding with value of 0...")
                        adata.obs[gene] = [np.nan] * adata.n_obs

        genes_up.sort()
        genes_down.sort()

        print(f"\tUsing upregulated genes {genes_up}...")
        print(f"\tUsing upregulated genes {genes_down}...")
        if adata.n_obs < 20:
            print("\tNot enough cells. Skipping cell type...")
        else:
            print("\tPlotting data...")
            fig = sc.pl.umap(
                adata,
                color=[
                    "leiden_res0_25",
                    "leiden_res0_50",
                    "leiden_res1_00",
                    "batch",
                    *genes_up,
                ],
                legend_loc="on data",
                show=False,
                return_fig=True,
            )
            fig.savefig(
                f"{output_folder_path}/umap_upregulated_{pathvalidate.sanitize_filename(sample)}_{pathvalidate.sanitize_filename(cell_type)}.svg"
            )
            fig = sc.pl.umap(
                adata,
                color=[
                    "leiden_res0_25",
                    "leiden_res0_50",
                    "leiden_res1_00",
                    "batch",
                    *genes_down,
                ],
                legend_loc="on data",
                show=False,
                return_fig=True,
            )
            fig.savefig(
                f"{output_folder_path}/umap_downregulated_{pathvalidate.sanitize_filename(sample)}_{pathvalidate.sanitize_filename(cell_type)}.svg"
            )
    else:
        print("\tFile does not exist. Skipping probably due to cell numbers being too low for clustering...")



def process_data(input_directory, output_folder_path):
    """
    Performs clustering.
    """
    print(f"Processing folder {input_directory}", flush=True)
    print("\tReading sample information...")
    info_path = os.path.join(input_directory, "info.csv")
    with open(info_path, newline="", encoding="utf-8") as csvfile:
        info_reader = csv.DictReader(
            csvfile, dialect="unix", delimiter=",", quotechar='"'
        )
        row = next(info_reader) # The file only has one row.
        sample_reference = row["reference sample"]
        sample_test = row["test sample"]
        cell_type = row["cell type"]

    print(f"\tProcessing samples {sample_reference} and {sample_test} for cell type {cell_type}")
    plot_genes(input_directory=input_directory, output_folder_path=output_folder_path, sample=sample_reference, cell_type=cell_type)
    plot_genes(input_directory=input_directory, output_folder_path=output_folder_path, sample=sample_test, cell_type=cell_type)


# Iterates over all sample directories and processes them conserving the directory structure.
for root, dirs, files in os.walk(INPUT_FOLDER_DGE):
    if "info.csv" in files:
        output_folder_path = os.path.join(
            MOUNT_PATHS["output"], root.removeprefix(INPUT_FOLDER_DGE)
        )
        os.makedirs(output_folder_path, exist_ok=True)
        process_data(root, output_folder_path)
