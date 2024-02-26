#!/usr/bin/python
"""This module calculates differentially expressed genes."""

import anndata
import anndata2ri
import csv
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

from rpy2.robjects.packages import importr

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER = MOUNT_PATHS["dependencies"]["clustered_pseudobulk_generation"] + "/"
FOLDER_COMPARISON_SAMPLE = "sample_comparison"
FOLDER_COMPARISON_CLUSTER = "cluster_comparison"
CLUSTER_KEY = "cluster"
SAMPLE_KEY = "sample"
REPLICATE_KEY = "replicate"
CLUSTERING_GROUP_KEY = "clusteringgroup"
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


def differential_gene_expression(
    pseudobulk_adata,
    sample_type_reference,
    sample_type_test,
    cluster,
    output_folder_path,
):
    """
    Performes differential gene expression analysis.
    """
    print(
        f"Comparing {sample_type_reference} and {sample_type_test} for cluster {cluster}",
        flush=True,
    )

    # Subsetting data.
    sample_reference_mask = pseudobulk_adata.obs[SAMPLE_KEY] == sample_type_reference
    sample_test_mask = pseudobulk_adata.obs[SAMPLE_KEY] == sample_type_test
    cluster_mask = pseudobulk_adata.obs[CLUSTER_KEY] == cluster
    data_mask = np.logical_and(
        np.logical_or(sample_reference_mask, sample_test_mask), cluster_mask
    )
    reference_cluster_mask = np.logical_and(sample_reference_mask, cluster_mask)
    test_cluster_mask = np.logical_and(sample_test_mask, cluster_mask)

    if pseudobulk_adata[reference_cluster_mask].n_obs < 2:
        print(
            f"\t{sample_type_reference} has only {pseudobulk_adata[reference_cluster_mask].n_obs} replicates. This is not enough to measure dispersion. Skipping cluster comparison..."
        )
        shutil.rmtree(output_folder_path)
        return

    if pseudobulk_adata[test_cluster_mask].n_obs < 2:
        print(
            f"\t{sample_type_test} has only {pseudobulk_adata[test_cluster_mask].n_obs} replicates. This is not enough to measure dispersion. Skipping cluster comparison..."
        )
        shutil.rmtree(output_folder_path)
        return

    adata_subset = pseudobulk_adata[data_mask]

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
        adata_subset, sample_type_reference, sample_type_test, output_folder_path
    )


def differential_gene_expression_exclude(
    pseudobulk_adata,
    pseudobulk_adata_exclude,
    clustering_group_test,
    cluster,
    output_folder_path,
):
    """
    Performes differential gene expression analysis.
    """
    print(
        f"Comparing {clustering_group_test} cluster {cluster} to the remaining sample...",
        flush=True,
    )

    # Subsetting data.
    sample_test_mask = pseudobulk_adata.obs[CLUSTERING_GROUP_KEY] == clustering_group_test
    sample_test_mask_exclude = (
        pseudobulk_adata_exclude.obs[CLUSTERING_GROUP_KEY] == clustering_group_test
    )

    cluster_mask = pseudobulk_adata.obs[CLUSTER_KEY] == cluster
    cluster_mask_exclude = pseudobulk_adata_exclude.obs[CLUSTER_KEY] == cluster

    data_mask = np.logical_and(sample_test_mask, cluster_mask)
    data_mask_exclude = np.logical_and(sample_test_mask_exclude, cluster_mask_exclude)

    adata_subset = pseudobulk_adata[data_mask]
    adata_subset_exclude = pseudobulk_adata_exclude[data_mask_exclude]

    if adata_subset.n_obs < 2:
        print(
            f"\t{clustering_group_test} has only {adata_subset.n_obs} replicates for cluster {cluster}. This is not enough to measure dispersion. Skipping cluster comparison..."
        )
        shutil.rmtree(output_folder_path)
        return

    if adata_subset_exclude.n_obs < 2:
        print(
            f"\t{clustering_group_test} has only {adata_subset_exclude.n_obs} replicates when excluding cluster {cluster}. This is not enough to measure dispersion. Skipping cluster comparison..."
        )
        shutil.rmtree(output_folder_path)
        return

    adata_merged = anndata.concat(
        [adata_subset, adata_subset_exclude],
        axis=0,
        join="outer",
        merge="same",
        uns_merge="same",
        label="excluded",
        keys=["no", "yes"],
    )

    print("\tLoading R dependencies...")
    importr("Seurat")
    importr("edgeR")
    print("\tRunning R...")
    edger_function = ro.r(
        """
        function(data, output_path){
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
            excluded <- colData(data)$excluded
            sample <- colData(data)$sample
            # Creates the design matrix.
            design <- model.matrix(~ 0 + excluded + sample)
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
            exclude_contrast <- makeContrasts("excludedno-excludedyes", levels = y$design)
            qlf <- glmQLFTest(fit, contrast=exclude_contrast)
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
        adata_merged, output_folder_path
    )


def get_adata(adata_file_name):
    """
    Reads the specified pseudobulk file,
    cleans the data and returns it.
    """
    adata_path = os.path.join(INPUT_FOLDER, adata_file_name)
    print(f"Reading pseudobulk data {adata_path}...")
    adata = anndata.read_h5ad(adata_path)

    # Replaces invalid characters in sample names.
    adata.obs[SAMPLE_KEY] = list(map(convert_string_to_r, adata.obs[SAMPLE_KEY]))
    adata.obs[SAMPLE_KEY] = adata.obs[SAMPLE_KEY].astype("category")

    # Replaces invalid characters in clustering group names.
    adata.obs[CLUSTERING_GROUP_KEY] = list(map(convert_string_to_r, adata.obs[CLUSTERING_GROUP_KEY]))
    adata.obs[CLUSTERING_GROUP_KEY] = adata.obs[CLUSTERING_GROUP_KEY].astype("category")

    return adata


# Parses the pseudobulk data.
adata_cluster = get_adata("cluster_pseudobulk.h5ad")
adata_cluster_exclude = get_adata("cluster_exclude_pseudobulk.h5ad")
adata_cluster_exclude.obs_names = list(map(lambda adata_obs_name: f"{adata_obs_name}_exclude",adata_cluster_exclude.obs_names))

print("Parsing information for sample comparisons...")
sample_comparison_path = os.path.join(MOUNT_PATHS["input"], "sample_comparison.csv")
sample_comparisons = []
with open(sample_comparison_path, newline="", encoding="utf-8") as csvfile:
    info_reader = csv.DictReader(csvfile, dialect="unix", delimiter=",", quotechar='"')
    for row in info_reader:
        sample_reference = row["reference sample"]
        sample_test = row["test sample"]
        sample_comparisons.append((sample_reference, sample_test))

# Runs differential gene expression analysis.
for sample_reference, sample_test in sample_comparisons:
    for cluster in adata_cluster.obs[CLUSTER_KEY].cat.categories:
        output_path = os.path.join(
            MOUNT_PATHS["output"],
            FOLDER_COMPARISON_SAMPLE,
            pathvalidate.sanitize_filename(f"{sample_reference}__vs__{sample_test}"),
            pathvalidate.sanitize_filename(cluster),
        )
        os.makedirs(output_path, exist_ok=True)
        with open(
            os.path.join(output_path, "info.csv"),
            mode="w",
            newline="",
            encoding="utf-8",
        ) as csvfile:
            info_writer = csv.writer(
                csvfile,
                dialect="unix",
                delimiter=",",
                quotechar='"',
                quoting=csv.QUOTE_MINIMAL,
            )
            info_writer.writerows(
                [
                    [
                        "reference sample",
                        "test sample",
                        "cluster",
                    ],
                    [
                        sample_reference,
                        sample_test,
                        cluster,
                    ],
                ]
            )
        differential_gene_expression(
            adata_cluster,
            sample_type_reference=convert_string_to_r(sample_reference),
            sample_type_test=convert_string_to_r(sample_test),
            cluster=cluster,
            output_folder_path=output_path,
        )


for clustering_group in adata_cluster.obs[CLUSTERING_GROUP_KEY].cat.categories:
    for cluster in adata_cluster.obs[CLUSTER_KEY].cat.categories:
        output_path = os.path.join(
            MOUNT_PATHS["output"],
            FOLDER_COMPARISON_CLUSTER,
            pathvalidate.sanitize_filename(clustering_group),
            pathvalidate.sanitize_filename(cluster),
        )
        os.makedirs(output_path, exist_ok=True)
        with open(
            os.path.join(output_path, "info.csv"),
            mode="w",
            newline="",
            encoding="utf-8",
        ) as csvfile:
            info_writer = csv.writer(
                csvfile,
                dialect="unix",
                delimiter=",",
                quotechar='"',
                quoting=csv.QUOTE_MINIMAL,
            )
            info_writer.writerows(
                [
                    [
                        "clustering_group",
                        "cluster",
                    ],
                    [
                        clustering_group,
                        cluster,
                    ],
                ]
            )
        differential_gene_expression_exclude(
            adata_cluster,
            adata_cluster_exclude,
            clustering_group_test=clustering_group,
            cluster=cluster,
            output_folder_path=output_path,
        )
