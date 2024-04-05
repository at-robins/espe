#!/usr/bin/python
"""This module runs the peak calling process."""

import json
import os
import subprocess

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER = next(iter(MOUNT_PATHS["dependencies"].values()))
GENOME_SIZE_PATH = os.path.join(
    MOUNT_PATHS["globals"]["GENOME"], "mappable_genome_size.txt"
)

with open(GENOME_SIZE_PATH, mode="rt", encoding="utf-8") as g_file:
    mappable_genome_size = g_file.read()
    print(f"Mappable genome size: {mappable_genome_size}")

options = (
    f"-f BAMPE -q 0.05 --keep-dup all --cutoff-analysis -B -g {mappable_genome_size}"
)

# Iterates over all sample directories and processes them conserving the directory structure.
FILE_EXTENSION = ".bam"
for root, dirs, files in os.walk(INPUT_FOLDER):
    for file in files:
        if file.endswith(FILE_EXTENSION):
            file_name = file.removesuffix(FILE_EXTENSION)
            input_file = os.path.join(root, file)
            base_output_path = os.path.join(
                MOUNT_PATHS["output"],
                root.removeprefix(INPUT_FOLDER + "/"),
            )
            print(f"Processing file {input_file}...", flush=True)
            os.makedirs(base_output_path, exist_ok=True)

            print("Calling peaks...", flush=True)
            subprocess.run(
                (
                    f"macs3 callpeak -t {input_file} "
                    f"-n {file_name}_narrow --outdir {base_output_path} "
                    f"{options}"
                ),
                cwd=base_output_path,
                shell=True,
                check=True,
            )
