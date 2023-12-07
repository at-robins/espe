#!/usr/bin/python
"""This module aggregates single cell data into pseudo-bulk data."""

import anndata
import json
import numpy as np
import os
import pandas as pd
import scanpy as sc

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER = MOUNT_PATHS["dependencies"]["ilc_composition"] + "/"
CELL_TYPE_KEY = "cell_type"
BATCH_KEY = "batch"
OUTPUT_CELL_TYPE_KEY = "celltype"
OUTPUT_SAMPLE_KEY = "sample"
OUTPUT_REPLICATE_KEY = "replicate"
OUTPUT_CELL_NUMBER = "cellcount"

# Setup of scanpy.
sc.settings.verbosity = 2


def aggregate_pseudobulk_data(file_path_sample, sample_name):
    """
    Aggregates pseudobulk data.
    """
    print(f"Processing file {file_path_sample}", flush=True)
    print("\tReading data...")
    adata = anndata.read_h5ad(file_path_sample)

    cell_types = adata.obs[CELL_TYPE_KEY].cat.categories
    batches = adata.obs[BATCH_KEY].cat.categories

    pseudobulk_dataframe = pd.DataFrame()
    cell_type_array = []
    sample_array = []
    replicate_array = []
    n_obs_array = []

    print(f"\tNumber of total observations: {adata.n_obs}")
    for cell_type in cell_types:
        for batch in batches:
            print(f"\tProcessing batch {batch} and cell type {cell_type}")
            cell_type_mask = adata.obs[CELL_TYPE_KEY] == cell_type
            batch_mask = adata.obs[BATCH_KEY] == batch
            adata_subset = adata[np.logical_and(cell_type_mask, batch_mask)]
            n_obs_subset = adata_subset.n_obs
            print(f"\t\tNumber of subset observations: {n_obs_subset}")
            if n_obs_subset > 0:
                aggregated_row = adata_subset.to_df().agg(np.sum)
                cell_type_array.append(cell_type)
                sample_array.append(sample_name)
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
        sample_array,
        replicate_array,
        cell_type_array,
        n_obs_array,
    )


# Aggregated pseudo-bulk scRNA-Seq data.
pseudobulk_data = pd.DataFrame()
samples = []
replicates = []
cell_types = []
n_cells = []

# Iterates over all sample directories and collects pseudo bulk data.
for root, dirs, files in os.walk(INPUT_FOLDER):
    for file in files:
        if file.casefold().endswith("filtered_feature_bc_matrix.h5ad"):
            input_file_path = os.path.join(root, file)
            (
                tmp_pseudobulk,
                tmp_samples,
                tmp_replicates,
                tmp_cell_types,
                tmp_n_cells,
            ) = aggregate_pseudobulk_data(
                file_path_sample=input_file_path,
                sample_name=root.removeprefix(INPUT_FOLDER),
            )
            pseudobulk_data = pd.concat(
                [pseudobulk_data, tmp_pseudobulk],
                ignore_index=False,
                join="outer",
            )
            samples.extend(tmp_samples)
            replicates.extend(tmp_replicates)
            cell_types.extend(tmp_cell_types)
            n_cells.extend(tmp_n_cells)

# Replaces NAs with zeros.
pseudobulk_data.fillna(0, inplace=True)

# Sets row names.
row_names = list(
    map(
        lambda name: name.replace(" ", "_"),
        [
            f"{sample}_{replicate}_{cell_type}"
            for sample, replicate, cell_type in zip(samples, replicates, cell_types)
        ],
    )
)
pseudobulk_data.index = row_names

# Sets additional variables.
pseudobulk_adata = anndata.AnnData(X=pseudobulk_data)
pseudobulk_adata.obs[OUTPUT_SAMPLE_KEY] = samples
pseudobulk_adata.obs[OUTPUT_REPLICATE_KEY] = replicates
pseudobulk_adata.obs[OUTPUT_CELL_TYPE_KEY] = cell_types
pseudobulk_adata.obs[OUTPUT_CELL_NUMBER] = n_cells

# Sets categorical data to a categorical data type.
pseudobulk_adata.obs[OUTPUT_SAMPLE_KEY] = pseudobulk_adata.obs[OUTPUT_SAMPLE_KEY].astype("category")
pseudobulk_adata.obs[OUTPUT_REPLICATE_KEY] = pseudobulk_adata.obs[OUTPUT_REPLICATE_KEY].astype("category")
pseudobulk_adata.obs[OUTPUT_CELL_TYPE_KEY] = pseudobulk_adata.obs[OUTPUT_CELL_TYPE_KEY].astype("category")

print("\tWriting pseudobulk data to file...")
pseudobulk_adata.write(f"{MOUNT_PATHS['output']}/pseudobulk.h5ad", compression="gzip")
