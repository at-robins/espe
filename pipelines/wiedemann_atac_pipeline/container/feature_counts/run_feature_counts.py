#!/usr/bin/python
"""This module creates the feature count matrix."""

import csv
import json
import math
import multiprocessing
import os
import subprocess
import sys

from pathlib import PurePath

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER_BAM = MOUNT_PATHS["dependencies"]["blacklist_removal"]
INPUT_FOLDER_PEAK = MOUNT_PATHS["dependencies"]["peak_annotation"]
BAM_SUFFIX = ".bam"
PEAK_SUFFIX = ".mergedPeak"
MATRIX_NON_SAMPLE_COLUMNS = 6

SAMPLE_INFO_DIRECTORY = "sample directory"
SAMPLE_INFO_TYPE = "sample type"
REPLICATE_KEY = "replicate_name"
SAMPLE_TYPE_KEY = "sample_type"

threads = math.floor(multiprocessing.cpu_count() * 0.8)
if threads < 1:
    threads = 1

# Sets batch information based on the supplied sample information.
print("Parsing sample information...")
sample_info_path = f"{MOUNT_PATHS['input']}/sample_information.csv"
sample_info_map = {}
with open(sample_info_path, newline="", encoding="utf-8") as csvfile:
    info_reader = csv.DictReader(csvfile, dialect="unix", delimiter=",", quotechar='"')
    for row in info_reader:
        sample_directory = row["sample directory"]
        sample_type = row["sample type"]
        sample_info_map[sample_directory] = sample_type

# Collects all peak files.
peak_files = []
for root, dirs, files in os.walk(INPUT_FOLDER_PEAK):
    for file in files:
        if file.endswith(PEAK_SUFFIX):
            peak_file_path = os.path.join(root, file)
            saf_file_path = os.path.join(MOUNT_PATHS["output"], "SAF", file)
            print(f"Detected peak file {peak_file_path}...", flush=True)
            print(f"Creating SAF file at {saf_file_path}...", flush=True)
            os.makedirs(os.path.dirname(saf_file_path), exist_ok=True)

            with open(
                peak_file_path, newline="", mode="rt", encoding="utf-8"
            ) as peak_file:
                with open(
                    saf_file_path, newline="", mode="wt", encoding="utf-8"
                ) as saf_file:
                    peak_reader = csv.reader(
                        peak_file, dialect="unix", delimiter="\t", quotechar='"'
                    )
                    saf_writer = csv.writer(
                        saf_file,
                        dialect="unix",
                        delimiter="\t",
                        quotechar='"',
                        quoting=csv.QUOTE_MINIMAL,
                    )
                    for i, row in enumerate(peak_reader):
                        if i == 0:
                            # Replaces the header.
                            saf_writer.writerow(
                                [
                                    "GeneID",
                                    "Chr",
                                    "Start",
                                    "End",
                                    "Strand",
                                ]
                            )
                        else:
                            # Copies relevant values of each other row.
                            saf_writer.writerow(
                                [
                                    row[0],
                                    row[1],
                                    row[2],
                                    row[3],
                                    row[4],
                                ]
                            )

            peak_files.append(saf_file_path)

if len(peak_files) != 1:
    print(
        f"One peak file was expected found {len(peak_files)}: {peak_files} ", flush=True
    )
    sys.exit(1)

# Collects all input BAM files.
bam_files = []
sample_dirs = []
for root, dirs, files in os.walk(INPUT_FOLDER_BAM):
    for file in files:
        if file.endswith(BAM_SUFFIX):
            bam_file_path = os.path.join(root, file)
            print(f"Detected BAM file {bam_file_path}...", flush=True)
            bam_files.append(bam_file_path)
            sample_dirs.append(root.removeprefix(INPUT_FOLDER_BAM).removeprefix("/"))


def rename_count_matrix(out_file):
    """
    Renames the sample headers of the count marix file.
    """
    in_file = f"{out_file}_tmp"
    os.rename(out_file, in_file)
    with open(in_file, newline="", mode="rt", encoding="utf-8") as count_matrix_in:
        with open(
            out_file, newline="", mode="wt", encoding="utf-8"
        ) as count_matrix_out:
            count_matrix_reader = csv.reader(
                count_matrix_in, dialect="unix", delimiter="\t", quotechar='"'
            )
            count_matrix_writer = csv.writer(
                count_matrix_out,
                dialect="unix",
                delimiter="\t",
                quotechar='"',
                quoting=csv.QUOTE_MINIMAL,
            )
            for i_matrix, row_matrix in enumerate(count_matrix_reader):
                # The first row is a comment, so we start with the second row.
                if i_matrix == 1 and len(row_matrix) > MATRIX_NON_SAMPLE_COLUMNS:
                    sample_rows = row_matrix[
                        MATRIX_NON_SAMPLE_COLUMNS : len(row_matrix)
                    ]
                    renamed_samples = []
                    for sample_path in sample_rows:
                        sample_folder = PurePath(sample_path).parent.name
                        sample_name = sample_info_map.get(sample_folder)
                        if sample_name is None:
                            print(
                                (
                                    "No sample information has been provided for directory "
                                    f"{sample_folder} in path {sample_path}."
                                ),
                                flush=True,
                            )
                            sys.exit(1)
                        renamed_samples.append(sample_name)
                    count_matrix_writer.writerow(
                        [*row_matrix[0:MATRIX_NON_SAMPLE_COLUMNS], *renamed_samples]
                    )
                else:
                    count_matrix_writer.writerow(row_matrix)
    os.remove(in_file)


print("Counting features...", flush=True)
for peak_file in peak_files:
    output_path = os.path.join(MOUNT_PATHS["output"], "counts", "count_matrix.txt")
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    subprocess.run(
        (
            f"featureCounts -F SAF -T {threads} -O -p --countReadPairs "
            f"-a {peak_file} -o {output_path} {' '.join(bam_files)}"
        ),
        shell=True,
        check=True,
    )
    rename_count_matrix(output_path)
