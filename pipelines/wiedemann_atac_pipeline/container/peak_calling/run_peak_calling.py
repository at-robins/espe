#!/usr/bin/python
"""This module runs the peak calling process."""

import json
import os
import subprocess

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER = MOUNT_PATHS["dependencies"]["blacklist_removal"]
GENOME_SIZE_PATH = MOUNT_PATHS["dependencies"]["mappable_size"]


# Iterates over all sample directories and processes them conserving the directory structure.
FILE_EXTENSION = ".bam"
for root, dirs, files in os.walk(INPUT_FOLDER):
    for file in files:
        if file.endswith(FILE_EXTENSION):
            relative_sub_folder = os.path.normpath(os.path.relpath(root, INPUT_FOLDER))
            file_name = file.removesuffix(FILE_EXTENSION)
            input_file = os.path.join(root, file)
            base_output_path = os.path.join(
                MOUNT_PATHS["output"],
                relative_sub_folder,
            )
            print(f"Processing file {input_file}...", flush=True)
            mappable_genome_size_file = os.path.join(
                GENOME_SIZE_PATH, relative_sub_folder, "unique_kmer_report.txt"
            )
            with open(mappable_genome_size_file, mode="rt", encoding="utf-8") as g_file:
                mappable_genome_size = g_file.readline().split()[0]
            print(f"\tMappable genome size: {mappable_genome_size}")
            options = f"-f BAMPE -q 0.05 --keep-dup all --cutoff-analysis -g {mappable_genome_size}"
            broad_output_folder = os.path.join(base_output_path, "broad")
            narrow_output_folder = os.path.join(base_output_path, "narrow")
            os.makedirs(broad_output_folder, exist_ok=True)
            os.makedirs(narrow_output_folder, exist_ok=True)

            print("\tCalling narrow peaks...", flush=True)
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

            print("\tCalling broad peaks...", flush=True)
            subprocess.run(
                (
                    f"macs3 callpeak -t {input_file} "
                    f"-n {file_name}_broad --outdir {broad_output_folder} "
                    f"--broad --broad-cutoff 0.1 "
                    f"{options}"
                ),
                cwd=broad_output_folder,
                shell=True,
                check=True,
            )
