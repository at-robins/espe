#!/usr/bin/python
"""This module runs the differential accessibility analysis."""

import csv
import json
import logging
import os
import pandas as pd
import pathvalidate
import rpy2.rinterface_lib.callbacks as rcb
import rpy2.robjects as ro

from pathlib import Path
from rpy2.robjects.packages import importr

# Setup of rpy2.
rcb.logger.setLevel(logging.INFO)
ro.r(
    '.libPaths(c("/usr/local/lib/R/site-library", "/usr/lib/R/site-library", "/usr/lib/R/library"))'
)

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER_COUNTS = MOUNT_PATHS["dependencies"]["feature_counts"]
INPUT_FOLDER_ANNOTATIONS = MOUNT_PATHS["dependencies"]["peak_annotation"]
COUNT_MATRIX_PATH = os.path.join(INPUT_FOLDER_COUNTS, "counts/count_matrix_final.txt")
ANNOTATION_PATH = os.path.join(INPUT_FOLDER_ANNOTATIONS, "merged.mergedPeak")
SAMPLE_COMPARISON_PATH = os.path.join(MOUNT_PATHS["input"], "sample_comparison.csv")
VALID_R_CHARACTERS = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_"


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


print("Loading count matrix...", flush=True)
peak_ids = []
count_matrix_values_by_row = []
groups = []
with open(
    COUNT_MATRIX_PATH, newline="", mode="rt", encoding="utf-8"
) as count_matrix_in:
    count_matrix_reader = csv.reader(
        count_matrix_in, dialect="unix", delimiter="\t", quotechar='"'
    )
    for i_matrix, row_matrix in enumerate(count_matrix_reader):
        peak_row = row_matrix[0]
        sample_rows = row_matrix[1 : len(row_matrix)]
        if i_matrix == 0:
            # Reads the condition names.
            groups = sample_rows
        else:
            # Reads the condition values.
            peak_ids.append(peak_row)
            for sample_value in sample_rows:
                count_matrix_values_by_row.append(int(sample_value))

print("Parsing information for sample comparisons...")
samples_test = []
samples_reference = []
sample_comparisons_output_suffixes = []
with open(SAMPLE_COMPARISON_PATH, newline="", encoding="utf-8") as csvfile:
    info_reader = csv.DictReader(csvfile, dialect="unix", delimiter=",", quotechar='"')
    for row in info_reader:
        sample_reference = row["reference sample"]
        sample_test = row["test sample"]
        samples_reference.append(sample_reference)
        samples_test.append(sample_test)
        sample_comparisons_output_suffixes.append(
            pathvalidate.sanitize_filename(f"{sample_test}__vs__{sample_reference}")
        )

r_count_matrix = ro.r.matrix(
    ro.IntVector(count_matrix_values_by_row),
    nrow=len(peak_ids),
    ncol=len(groups),
    byrow=True,
)

print("Loading R dependencies...", flush=True)
importr("ggplot2")
importr("DESeq2")
print("Running R...", flush=True)
deseq_function = ro.r(
    """
    function(
        data_matrix,
        peak_ids,
        groups,
        group_names,
        conditions_test,
        conditions_reference,
        condition_names_test,
        condition_names_reference,
        output_path,
        output_suffixes
    ) {
        # Creates a mapping between groups and their names.
        names(group_names) <- groups
        cat("Preprocessing data matrix...\\n", sep="")
        # Names rows and columns of the data matrix.
        rownames(data_matrix) <- peak_ids
        colnames(data_matrix) <- 1:ncol(data_matrix)
        meta_data <- data.frame(
            factor(groups),
            row.names = colnames(data_matrix)
        )
        colnames(meta_data) <- c("groups")
        atac_dds <- DESeqDataSetFromMatrix(
            countData = data_matrix,
            colData = meta_data,
            design = ~groups
        )
        atac_dds <- DESeq(atac_dds, parallel = TRUE)

        cat("Regularising count matrix...\\n", sep="")
        atac_rlog <- rlog(atac_dds)

        cat("Exporting count matrix...\\n", sep="")
        atac_matrix_regularised <- as.data.frame(assay(atac_rlog))
        colnames(atac_matrix_regularised) = group_names
        write.csv(
            atac_matrix_regularised,
            file = paste(output_path, "count_matrix_regularised.csv", sep = "/"),
            row.names=TRUE
        )

        cat("Plotting PCA...\\n", sep="")
        pca_plot <- plotPCA(atac_rlog, intgroup = c("groups"), ntop = nrow(atac_rlog)) +
            scale_colour_discrete(labels = group_names, name = "Condition") +
            theme(
                panel.background = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line = element_blank(),
                panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.0)
            )
        ggsave(filename = paste(output_path, "pca.svg", sep = "/"), plot = pca_plot)

        cat("Performing differential accessiblity analysis...\\n", sep = "")
        for (i in 1:length(conditions_test)) {
            cat(
                "\\tComparing sample ",
                condition_names_test[i],
                " to reference ",
                condition_names_reference[i],
                "...\\n",
                sep=""
            )
            atac_result <- results(
                atac_dds,
                contrast = c("groups", conditions_test[i], conditions_reference[i])
            )
            atac_result["conditionSample"] <- condition_names_test[i]
            atac_result["conditionReference"] <- condition_names_reference[i]
            write.csv(
                atac_result,
                file = paste(
                    output_path,
                    "/",
                    "differential_accessibility_",
                    output_suffixes[i],
                    ".csv",
                    sep = ""
                ),
                row.names=TRUE
            )

            cat("Plotting MA...\\n", sep="")
            svg(
                file = paste(
                    output_path,
                    "/",
                    "ma_",
                    output_suffixes[i],
                    ".svg",
                    sep = ""
                )
            )
            plotMA(atac_result)
            dev.off()
        }
    }
    """
)
deseq_function(
    r_count_matrix,
    ro.StrVector(peak_ids),
    ro.StrVector(list(map(convert_string_to_r, groups))),
    ro.StrVector(groups),
    ro.StrVector(list(map(convert_string_to_r, samples_test))),
    ro.StrVector(list(map(convert_string_to_r, samples_reference))),
    ro.StrVector(samples_test),
    ro.StrVector(samples_reference),
    MOUNT_PATHS["output"],
    ro.StrVector(sample_comparisons_output_suffixes),
)

print("Loading annotations...", flush=True)
annotation = pd.read_csv(
    ANNOTATION_PATH, sep="\t", header=0, index_col=0, encoding="utf-8"
)

# Annotates all output files.
for root, dirs, files in os.walk(MOUNT_PATHS["output"]):
    for file in files:
        if "differential_accessibility" in file and file.endswith(".csv"):
            raw_path = Path(os.path.join(root, file))
            annotated_path = Path(os.path.join(root, f"annotated_{file}"))
            print(f"Annotating {raw_path}...", flush=True)
            da_table = pd.read_csv(
                raw_path, sep=",", header=0, index_col=0, encoding="utf-8"
            )
            pd.merge(
                left=da_table,
                right=annotation,
                left_index=True,
                right_index=True,
                how="left",
            ).sort_values("padj", ascending=True).to_csv(
                annotated_path,
                sep=",",
                encoding="utf-8",
                index=True,
                index_label="PeakID",
            )
