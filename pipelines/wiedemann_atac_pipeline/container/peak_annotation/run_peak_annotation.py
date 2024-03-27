#!/usr/bin/python
"""This module runs the peak annotation."""

import json
import os
import subprocess

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER = next(iter(MOUNT_PATHS["dependencies"].values()))
ANNOTATION_PATH = os.path.join(MOUNT_PATHS["globals"]["GENOME"], "annotations.gtf")

# Iterates over all sample directories and processes them conserving the directory structure.
print("Annotating peaks...", flush=True)
for root, dirs, files in os.walk(INPUT_FOLDER):
    for file in files:
        if file.startswith("combined_sorted."):
            file_base_name_output = file.replace(".", "_")
            file_path_input = os.path.join(root, file)
            directory_path_output = os.path.join(
                MOUNT_PATHS["output"],
                os.path.normpath(os.path.relpath(root, INPUT_FOLDER)),
            )
            print(f"\tProcessing file {file_path_input}...", flush=True)
            os.makedirs(directory_path_output, exist_ok=True)

            subprocess.run(
                (
                    f"uropa --bed {file_path_input} --gtf {ANNOTATION_PATH} "
                    f"-p {os.path.join(directory_path_output, file_base_name_output)} -o {directory_path_output} "
                    "--feature_anchor start --distance 20000 10000 --feature gene --show_attributes gene_name"
                ),
                cwd=directory_path_output,
                shell=True,
                check=True,
            )

print("Extracting headers...", flush=True)
for root, dirs, files in os.walk(MOUNT_PATHS["output"]):
    for file in files:
        if file.endswith("_finalhits.txt"):
            hit_file_path_input = os.path.join(root, file)
            hit_file_path_output = os.path.join(
                root, f"{file.removesuffix('.txt')}_header.txt"
            )
            print(f"\tExtracting header from file {file_path_input}...", flush=True)
            with open(hit_file_path_input, mode="rt", encoding="utf-8") as hit_file_in:
                with open(
                    hit_file_path_output, mode="wt", encoding="utf-8"
                ) as hit_file_out:
                    hit_file_out.write(next(hit_file_in))
