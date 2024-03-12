#!/usr/bin/python
"""This module runs the peak calling process."""

import json
import os
import sys

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER = MOUNT_PATHS["dependencies"]["blacklist_removal"]
CHROM_SIZES_PATH = os.path.join(MOUNT_PATHS["globals"]["GENOME"], "genome.chrom.sizes")

options = "-f BEDPE -q 0.01 -B -g `awk '{SUM+=$2}END{print SUM}' {CHROM_SIZES_PATH}`"

# Iterates over all sample directories and processes them conserving the directory structure.
FILE_EXTENSION = ".bedGraph"
INPUT_SUFFIX_FORWARD = f"_nucleosomefree{FILE_EXTENSION}"
for root, dirs, files in os.walk(INPUT_FOLDER):
    for file in files:
        if file.endswith(INPUT_SUFFIX_FORWARD):
            file_name = file.removesuffix(FILE_EXTENSION)
            input_file = os.path.join(root, file)
            base_output_path = os.path.join(
                MOUNT_PATHS["output"],
                root.removeprefix(INPUT_FOLDER + "/"),
            )
            print(f"Processing file {input_file}...", flush=True)
            broad_output_folder = os.path.join(base_output_path, "broad")
            narrow_output_folder = os.path.join(base_output_path, "narrow")
            os.makedirs(broad_output_folder, exist_ok=True)
            os.makedirs(narrow_output_folder, exist_ok=True)

            print("Calling narrow peaks...", flush=True)
            subprocess.run(
                (
                    f"macs3 callpeak -t {input_file} "
                    f"-n {file_name}_narrow --outdir {narrow_output_folder} "
                    f"{options}"
                ),
                cwd=narrow_output_folder,
                shell=True,
                check=True,
            )

            print("Calling broad peaks...", flush=True)
            subprocess.run(
                (
                    f"macs3 callpeak -t {input_file} "
                    f"-n {file_name}_broad --outdir {broad_output_folder} "
                    f"--broad --broad-cutoff 0.01 "
                    f"{options}"
                ),
                cwd=narrow_output_folder,
                shell=True,
                check=True,
            )
