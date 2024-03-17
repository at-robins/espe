#!/usr/bin/python
"""This module calculates cluster relation trees."""

import json
import os
import subprocess


MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER = MOUNT_PATHS["dependencies"]["cluster_resolution_data"] + "/"


def process_data(file_path_input, output_folder_path):
    """
    Computes cluster statistics.
    """
    print(f"Processing file {file_path_input}", flush=True)
    subprocess.run(
        f"leiden-optimisation -s 0.85 -o {output_folder_path} {file_path_input}",
        shell=True,
        check=True,
    )

print("Searching for data...")
# Iterates over all sample directories and processes them conserving the directory structure.
for root, dirs, files in os.walk(INPUT_FOLDER):
    for file in files:
        if file.casefold().endswith("cluster_resolution_data.csv"):
            file_path_input = os.path.join(root, file)
            output_folder_path = os.path.join(
                MOUNT_PATHS["output"], root.removeprefix(INPUT_FOLDER)
            )
            os.makedirs(output_folder_path, exist_ok=True)
            process_data(file_path_input, output_folder_path)
