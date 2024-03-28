#!/usr/bin/python
"""This module creates the feature count matrix."""

import csv
import json
import math
import multiprocessing
import os
import subprocess

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER_BAM = MOUNT_PATHS["dependencies"]["blacklist_removal"]
INPUT_FOLDER_PEAK = MOUNT_PATHS["dependencies"]["peak_annotation"]
BAM_SUFFIX = ".bam"
PEAK_NARROW_SUFFIX = ".narrowPeak"
PEAK_BROAD_SUFFIX = ".broadPeak"


threads = math.floor(multiprocessing.cpu_count() * 0.8)
if threads < 1:
    threads = 1

peak_files = []
for root, dirs, files in os.walk(INPUT_FOLDER_PEAK):
    for file in files:
        if file.endswith(PEAK_BROAD_SUFFIX) or file.endswith(PEAK_NARROW_SUFFIX):
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
                            saf_writer.writerow(
                                [
                                    "GeneID",
                                    row[1],
                                    row[2],
                                    row[3],
                                    row[4],
                                ]
                            )
                        else:
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

print("Counting features...", flush=True)
# Iterates over all sample directories and processes them conserving the directory structure.
for root, dirs, files in os.walk(INPUT_FOLDER_BAM):
    for file in files:
        if file.endswith(f"nucleosomefree{BAM_SUFFIX}"):
            print(f"\tProcessing file {os.path.join(root, file)}...", flush=True)
            file_input_path = os.path.join(root, file)
            file_base_name = file.removesuffix(BAM_SUFFIX)
            file_base_path = os.path.join(root, file_base_name)
            file_base_output_path = os.path.join(
                MOUNT_PATHS["output"],
                file_base_path.removeprefix(INPUT_FOLDER_BAM + "/"),
            )
            os.makedirs(os.path.dirname(file_base_output_path), exist_ok=True)
            for peak_file in peak_files:
                if peak_file.endswith(PEAK_NARROW_SUFFIX):
                    file_output_path = f"{file_base_output_path}_narrow.txt"
                else:
                    file_output_path = f"{file_base_output_path}_broad.txt"

                subprocess.run(
                    (
                        f"featureCounts -F SAF -T {threads} -O -p --countReadPairs "
                        f"-a {peak_file} -o {file_output_path} {file_input_path}"
                    ),
                    shell=True,
                    check=True,
                )
