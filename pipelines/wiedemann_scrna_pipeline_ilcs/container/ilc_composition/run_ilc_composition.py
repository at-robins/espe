#!/usr/bin/python
"""This module performs analyses the group 1 ILC composition."""

import anndata
import csv
import json
import logging
import numpy as np
import os
import scanpy as sc
import seaborn as sns

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER = MOUNT_PATHS["input"] + "/"

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


def process_data(file_path_input, output_folder_path, metrics_writer):
    """
    Performs analyses.
    """
    print(f"Processing file {file_path_input}", flush=True)
    print("\tReading data...")
    adata = anndata.read_h5ad(file_path_input)
    # adata.layers["non_normalised"] = adata.X

    # layer_used = "non_normalised"
    layer_used = "log1p_norm"
    # layer_used = "scran_normalisation"
    # layer_used = "analytic_pearson_residuals"

    n_genes_plot_raw = sc.pl.scatter(
        adata=adata,
        x="Eomes",
        y="Zfp683",
        layers=layer_used,
        title="EOMES vs. HOBIT",
        show=False,
    )
    n_genes_plot_raw.set(xlabel="EOMES", ylabel="HOBIT")
    n_genes_plot_raw.figure.savefig(f"{output_folder_path}/scatter_eomes_hobit.svg")

    eomes = adata[:, "Eomes"].layers[layer_used].toarray()
    hobit = adata[:, "Zfp683"].layers[layer_used].toarray()
    cd49a = adata[:, "Itga1"].layers[layer_used].toarray()
    cd49b = adata[:, "Itga2"].layers[layer_used].toarray()
    tbet = adata[:, "Tbx21"].layers[layer_used].toarray()
    nkp46 = adata[:, "Ncr1"].layers[layer_used].toarray()
    nk11 = adata[:, "Klrb1"].layers[layer_used].toarray()

    eomes_positive_mask = np.logical_and(eomes > 0, hobit <= 0)
    hobit_positive_mask = np.logical_and(hobit > 0, eomes <= 0)
    double_positive_mask = np.logical_and(eomes > 0, hobit > 0)
    double_negative_mask = np.logical_and(eomes <= 0, hobit <= 0)

    adata.obs["cell_type"] = np.select(
        [
            eomes_positive_mask,
            hobit_positive_mask,
            double_positive_mask,
            double_negative_mask,
        ],
        ["NK cell", "ILC1", "intermediate group 1 ILC", "negative"],
    )
    # sc.tl.rank_genes_groups(
    #     adata, groupby="cell_type", layer=layer_used, method="wilcoxon"
    # )
    # print(
    #     sc.get.rank_genes_groups_df(adata, group="negative").to_string(),
    #     flush=True,
    # )

    print("\tWriting metrics to file...")
    metrics_writer.writerow(
        [
            file_path_input.removeprefix(INPUT_FOLDER),
            np.count_nonzero(eomes_positive_mask),
            np.count_nonzero(hobit_positive_mask),
            np.count_nonzero(double_positive_mask),
            np.count_nonzero(double_negative_mask),
        ]
    )

    print("\tWriting data to file...")
    adata.write(f"{output_folder_path}/filtered_feature_bc_matrix.h5ad", compression="gzip")

    print("\tPlotting data...")
    sns.displot(data=eomes, bins=100, kde=False).set(
        xlabel="EOMES expression", ylabel="Number of cells"
    ).savefig(f"{output_folder_path}/histo_eomes.svg")

    sns.displot(data=hobit, bins=100, kde=False).set(
        xlabel="HOBIT expression", ylabel="Number of cells"
    ).savefig(f"{output_folder_path}/histo_hobit.svg")

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


with open(
    f"{MOUNT_PATHS['output']}/metrics.csv", mode="w", newline="", encoding="utf-8"
) as csvfile:
    metrics_writer = csv.writer(
        csvfile, dialect="unix", delimiter=",", quotechar='"', quoting=csv.QUOTE_MINIMAL
    )
    metrics_writer.writerow(
        [
            "Sample",
            "Total number of cells",
            "Number of NK cells",
            "Number of ILC1s",
            "Number of intermediate group 1 ILCs",
            "Number of negative cells",
        ]
    )
    # Iterates over all sample directories and processes them conserving the directory structure.
    for root, dirs, files in os.walk(INPUT_FOLDER):
        for file in files:
            if file.casefold().endswith("filtered_feature_bc_matrix.h5ad"):
                file_path_input = os.path.join(root, file)
                output_folder_path = os.path.join(
                    MOUNT_PATHS["output"], root.removeprefix(INPUT_FOLDER)
                )
                os.makedirs(output_folder_path, exist_ok=True)
                process_data(file_path_input, output_folder_path, metrics_writer)
