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
import scanpy as sc
import seaborn as sns
import rpy2.rinterface_lib.callbacks as rcb
import rpy2.robjects as ro
import warnings

from rpy2.robjects.packages import importr

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER = MOUNT_PATHS["dependencies"]["pseudobulk_generation"] + "/"
CELL_TYPE_KEY = "celltype"
SAMPLE_KEY = "sample"
REPLICATE_KEY = "replicate"

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
    return val.replace(" ", "_")


def differential_gene_expression(
    pseudobulk_adata,
    sample_type_reference,
    sample_type_test,
    cell_type,
    output_folder_path,
):
    """
    Performes differential gene expression analysis.
    """
    print(
        f"Comparing {sample_type_reference} and {sample_type_test} for cell type {cell_type}",
        flush=True,
    )

    # Subsetting data.
    sample_reference_mask = adata.obs[SAMPLE_KEY] == sample_type_reference
    sample_test_mask = adata.obs[SAMPLE_KEY] == sample_type_test
    cell_type_mask = adata.obs[CELL_TYPE_KEY] == cell_type
    data_mask = np.logical_and(
        np.logical_or(sample_reference_mask, sample_test_mask), cell_type_mask
    )
    adata_subset = adata[data_mask]

    # Performs a preliminary PCA.
    # adata_subset.layers['counts'] = adata_subset.X.copy()
    # sc.pp.normalize_total(adata_subset, target_sum=1e6)
    # sc.pp.log1p(adata_subset)
    # sc.pp.pca(adata_subset)
    # adata_subset.obs["lib_size"] = np.sum(adata_subset.layers["counts"], axis=1)
    # adata_subset.obs["log_lib_size"] = np.log(adata_subset.obs["lib_size"])
    # sc.pl.pca(adata_subset, color=adata_subset.obs, ncols=1, size=300, save=True)
    # adata_subset.X = adata_subset.layers['counts'].copy()

    print("\tLoading R dependencies...")
    importr("Seurat")
    importr("edgeR")
    print("\tRunning R...")
    edger_function = ro.r(
        """
        function(data, sample_reference, sample_test, output_path){
            # create an edgeR object with counts and grouping factor
            y <- DGEList(assay(data, "X"), group = colnames(data))
            # filter out genes with low counts
            print("Dimensions before subsetting:")
            print(dim(y))
            keep <- filterByExpr(y)
            y <- y[keep, , keep.lib.sizes=FALSE]
            print("Dimensions after subsetting:")
            print(dim(y))
            # normalize
            y <- calcNormFactors(y)
            # create a vector that is concatentation of condition and cell type that we will later use with contrasts
            sample <- colData(data)$sample
            # create a design matrix: here we have multiple donors so also consider that in the design matrix
            design <- model.matrix(~ 0 + sample)
            # estimate dispersion
            y <- estimateDisp(y, design = design)
            # fit the model
            fit <- glmQLFit(y, design)
            print("Plotting data...")
            svg(paste(output_path, "mds_plot.svg", sep = "/"))
            plotMDS(y)
            dev.off()
            svg(paste(output_path, "bcv_plot.svg", sep = "/"))
            plotBCV(y)
            dev.off()
            # eval workaround as make makeContrasts does not accept a string variable.
            contrast_string = paste("sample", sample_test, " - sample", sample_reference, sep = "")
            cmd <- paste("myContrast <- makeContrasts(", contrast_string, ", levels = y$design)", sep ='"')
            eval(parse(text = cmd))
            qlf <- glmQLFTest(fit, contrast=myContrast)
            # get all of the DE genes and calculate Benjamini-Hochberg adjusted FDR
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


print("Reading pseudobulk data...")
adata = anndata.read_h5ad(f"{INPUT_FOLDER}pseudobulk.h5ad")

# Replaces invalid characters in sample names.
adata.obs[SAMPLE_KEY] = list(map(convert_string_to_r, adata.obs[SAMPLE_KEY]))
adata.obs[SAMPLE_KEY] = adata.obs[SAMPLE_KEY].astype("category")


print("Parsing information for sample comparisons...")
sample_comparison_path = f"{MOUNT_PATHS['input']}/sample_comparison.csv"
sample_comparisons = []
with open(sample_comparison_path, newline="", encoding="utf-8") as csvfile:
    info_reader = csv.DictReader(csvfile, dialect="unix", delimiter=",", quotechar='"')
    for row in info_reader:
        sample_reference = row["reference sample"]
        sample_test = row["test sample"]
        sample_comparisons.append((sample_reference, sample_test))

# Runs differential gene expression analysis.
for sample_reference, sample_test in sample_comparisons:
    for cell_type in adata.obs[CELL_TYPE_KEY].cat.categories:
        output_path = os.path.join(
            MOUNT_PATHS["output"],
            pathvalidate.sanitize_filename(f"{sample_reference}__vs__{sample_test}"),
            pathvalidate.sanitize_filename(cell_type),
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
                        "cell type",
                    ],
                    [
                        sample_reference,
                        sample_test,
                        cell_type,
                    ],
                ]
            )
        differential_gene_expression(
            adata,
            sample_type_reference=convert_string_to_r(sample_reference),
            sample_type_test=convert_string_to_r(sample_test),
            cell_type=cell_type,
            output_folder_path=output_path,
        )
