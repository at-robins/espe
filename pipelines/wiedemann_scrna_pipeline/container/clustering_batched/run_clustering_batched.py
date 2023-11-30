#!/usr/bin/python
"""This module performs data clustering for QC."""

import anndata
import csv
import json
import logging
import numpy as np
import os
import scanpy as sc
import seaborn as sns

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER = MOUNT_PATHS["dependencies"]["feature_selection"] + "/"
UNBATCHED_SUBFOLDER = "batched"
INPUT_FOLDER_BATCHED = f"{INPUT_FOLDER}{UNBATCHED_SUBFOLDER}/"

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


def process_data(file_path_input, output_folder_path):
    """
    Performs clustering.
    """
    print(f"Processing file {file_path_input}", flush=True)
    print("\tReading filtered data...")
    adata = anndata.read_h5ad(file_path_input)

    # print("\tPerforming clustering...")
    # adata.X = adata.layers["log1p_norm"]
    # # Performing HGV and PCA first to reduce dimensionality for UMAP.
    # adata.var["highly_variable"] = adata.var["highly_deviant"]
    # sc.pp.pca(adata, svd_solver="arpack", use_highly_variable=True)
    # sc.pp.neighbors(adata, n_pcs=30)
    # sc.tl.umap(adata)
    # sc.tl.leiden(adata, key_added="leiden_res0_25", resolution=0.25)
    # sc.tl.leiden(adata, key_added="leiden_res0_5", resolution=0.5)
    # sc.tl.leiden(adata, key_added="leiden_res1", resolution=1.0)

    # print("\tPlotting data...")
    # fig = sc.pl.umap(
    #     adata,
    #     color=["leiden_res0_25", "leiden_res0_5", "leiden_res1"],
    #     legend_loc="on data",
    #     show=False,
    #     return_fig=True,
    # )
    # fig.tight_layout()
    # fig.savefig(f"{output_folder_path}/umap.svg")

    n_genes_plot_raw = sc.pl.scatter(
        adata=adata,
        x="Eomes",
        y="Zfp683",
        layers = "log1p_norm",
        title="EOMES vs HOBIT",
        show=False,
    )
    n_genes_plot_raw.set(
        xlabel="EOMES", ylabel="HOBIT"
    )
    n_genes_plot_raw.figure.savefig(
        f"{output_folder_path}/scatter_eomes_hobit.svg"
    )

    eomes = adata[:, "Eomes"].layers["log1p_norm"].toarray()
    hobit = adata[:, "Zfp683"].layers["log1p_norm"].toarray()
    cd49a = adata[:, "Itga1"].layers["log1p_norm"].toarray()
    cd49b = adata[:, "Itga2"].layers["log1p_norm"].toarray()
    tbet = adata[:, "Tbx21"].layers["log1p_norm"].toarray()
    nkp46 = adata[:, "Ncr1"].layers["log1p_norm"].toarray()
    nk11 = adata[:, "Klrb1"].layers["log1p_norm"].toarray()

    eomes_positve_mask = np.logical_and(eomes > 0, hobit <= 0)
    hobit_positive_mask = np.logical_and(hobit > 0, eomes <= 0)
    double_positive_mask = np.logical_and(eomes > 0, hobit > 0)
    double_negative_mask = np.logical_and(eomes <= 0, hobit <= 0)

    adata.obs["double_negative"] = np.select([double_negative_mask, ~double_negative_mask], ["double negative", "positive"])
    sc.tl.rank_genes_groups(adata, groupby = "double_negative", layer="log1p_norm", method='wilcoxon')
    print(sc.get.rank_genes_groups_df(adata, group = "double negative").to_string(), flush=True)

    print(f"EOMES positive: {len(eomes[np.logical_and(eomes > 0, hobit == 0)]) / len(eomes)}")
    print(f"HOBIT positive: {len(hobit[np.logical_and(hobit > 0, eomes == 0)]) / len(eomes)}")
    print(f"double positive: {len(eomes[np.logical_and(eomes > 0, hobit > 0)]) / len(eomes)}")
    print(f"double negative: {len(eomes[np.logical_and(eomes == 0, hobit == 0)]) / len(eomes)}")

    sns.displot(data=eomes, bins=100, kde=False).set(
        xlabel="EOMES expression", ylabel="Number of cells"
    ).savefig(f"{output_folder_path}/histo_eomes.svg")

    sns.displot(data=eomes[hobit > 0], bins=100, kde=False).set(
        xlabel="EOMES expression", ylabel="Number of cells"
    ).savefig(f"{output_folder_path}/histo_eomes_no_0hobit.svg")

    sns.displot(data=hobit, bins=100, kde=False).set(
        xlabel="HOBIT expression", ylabel="Number of cells"
    ).savefig(f"{output_folder_path}/histo_hobit.svg")

    sns.displot(data=hobit[eomes > 0], bins=100, kde=False).set(
        xlabel="HOBIT expression", ylabel="Number of cells"
    ).savefig(f"{output_folder_path}/histo_hobit_no_0eomes.svg")

    sns.displot(data=cd49a[double_negative_mask], bins=100, kde=False).set(
        xlabel="CD49A expression", ylabel="Number of cells"
    ).savefig(f"{output_folder_path}/histo_cd49a_double_negative.svg")

    sns.displot(data=cd49b[double_negative_mask], bins=100, kde=False).set(
        xlabel="CD49B expression", ylabel="Number of cells"
    ).savefig(f"{output_folder_path}/histo_cd49b_double_negative.svg")

    sns.displot(data=tbet[double_negative_mask], bins=100, kde=False).set(
        xlabel="TBET expression", ylabel="Number of cells"
    ).savefig(f"{output_folder_path}/histo_tbet_double_negative.svg")

    sns.displot(data=nk11[double_negative_mask], bins=100, kde=False).set(
        xlabel="NK1.1 expression", ylabel="Number of cells"
    ).savefig(f"{output_folder_path}/histo_nk11_double_negative.svg")

    sns.displot(data=nkp46[double_negative_mask], bins=100, kde=False).set(
        xlabel="NKP46 expression", ylabel="Number of cells"
    ).savefig(f"{output_folder_path}/histo_nkp46_double_negative.svg")


# Iterates over all sample directories and processes them conserving the directory structure.
for root, dirs, files in os.walk(INPUT_FOLDER_BATCHED):
    for file in files:
        if file.casefold().endswith("filtered_feature_bc_matrix.h5ad"):
            file_path_input = os.path.join(root, file)
            output_folder_path = os.path.join(
                MOUNT_PATHS["output"], root.removeprefix(INPUT_FOLDER_BATCHED)
            )
            os.makedirs(output_folder_path, exist_ok=True)
            process_data(file_path_input, output_folder_path)
