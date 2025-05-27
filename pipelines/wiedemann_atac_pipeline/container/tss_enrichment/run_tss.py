#!/usr/bin/python
"""This module plots transcription start site enrichment."""

import json
import math
import multiprocessing
import os
import subprocess
import sys


MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER = next(iter(MOUNT_PATHS["dependencies"].values()))
FILE_EXTENSION_BAM = ".bam"
FILE_EXTENSION_BIGWIG = ".bw"
FILE_EXTENSION_MATRIX = ".mat.gz"
FILE_EXTENSION_PLOT = ".svg"
GENE_ANNOTATION_PATH = os.path.join(MOUNT_PATHS["globals"]["GENOME"], "annotations.gtf")


threads = math.floor(multiprocessing.cpu_count() * 0.8)
if threads < 1:
    threads = 1

print("Generating bigWig files...", flush=True)
for root, dirs, files in os.walk(INPUT_FOLDER):
    for file in files:
        if file.endswith(FILE_EXTENSION_BAM):
            bam_file = os.path.join(root, file)
            print(f"\tProcessing {bam_file}...", flush=True)

            output_path_folder = os.path.join(
                MOUNT_PATHS["output"],
                os.path.normpath(os.path.relpath(root, INPUT_FOLDER)),
            )

            output_path_file = os.path.join(
                output_path_folder,
                file.removesuffix(FILE_EXTENSION_BAM) + FILE_EXTENSION_BIGWIG,
            )

            os.makedirs(output_path_folder, exist_ok=True)

            print("\tConverting to bigWig...", flush=True)
            subprocess.run(
                (
                    "bamCoverage "
                    f'--bam "{bam_file}" '
                    f'--outFileName "{output_path_file}" '
                    f"--numberOfProcessors {threads} "
                    "--outFileFormat bigwig"
                ),
                cwd=output_path_folder,
                shell=True,
                check=True,
            )

print("Generating transcription site enrichment plots...", flush=True)
for root, dirs, files in os.walk(MOUNT_PATHS["output"]):
    for file in files:
        if file.endswith(FILE_EXTENSION_BIGWIG):
            bw_file = os.path.join(root, file)
            print(f"\tProcessing {bw_file}...", flush=True)

            output_path_folder = root

            print("\tCalculating count matrix...", flush=True)
            output_path_matrix = os.path.join(
                output_path_folder,
                file.removesuffix(FILE_EXTENSION_BIGWIG) + FILE_EXTENSION_MATRIX,
            )

            subprocess.run(
                (
                    "computeMatrix reference-point "
                    "--referencePoint TSS "
                    "--afterRegionStartLength 1000 "
                    "--beforeRegionStartLength 1000 "
                    "--skipZeros "
                    "--exonID exon "
                    "--transcriptID transcript "  # TODO: fix: transcript
                    f"--numberOfProcessors {threads} "
                    f'--regionsFileName "{GENE_ANNOTATION_PATH}" '
                    f'--scoreFileName "{bw_file}" '
                    f'--outFileName "{output_path_matrix}"'
                ),
                cwd=output_path_folder,
                shell=True,
                check=True,
            )

            print("\tPlotting read density...", flush=True)
            output_path_line_plot = os.path.join(
                output_path_folder,
                bw_file.removesuffix(FILE_EXTENSION_BIGWIG)
                + "_read_density"
                + FILE_EXTENSION_PLOT,
            )

            subprocess.run(
                (
                    "plotProfile "
                    "--refPointLabel TSS "
                    '--plotTitle "Read density" '
                    f"-m {output_path_matrix} "
                    f'-o "{output_path_line_plot}" '
                ),
                cwd=output_path_folder,
                shell=True,
                check=True,
            )

            print("\tPlotting heatmap...", flush=True)
            output_path_heatmap = os.path.join(
                output_path_folder,
                bw_file.removesuffix(FILE_EXTENSION_BIGWIG)
                + "_heatmap"
                + FILE_EXTENSION_PLOT,
            )

            subprocess.run(
                (
                    "plotHeatmap "
                    "--refPointLabel TSS "
                    '--plotTitle "Read distribution heatmap" '
                    '--whatToShow "heatmap and colorbar" '
                    f"--matrixFile {output_path_matrix} "
                    f'--outFileName "{output_path_heatmap}" '
                ),
                cwd=output_path_folder,
                shell=True,
                check=True,
            )
