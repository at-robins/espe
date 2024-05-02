#!/usr/bin/python
"""This module shows marker gene expression in the clustered data."""

import json
import os
import subprocess

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER = next(iter(MOUNT_PATHS["dependencies"].values()))


def process_data(file_path_input, output_folder_path):
    """
    .
    """
    print(f"Processing file {file_path_input}", flush=True)
    subprocess.run(
        (
            "python /TFvelo/TFvelo_run_demo.py "
            f"--dataset_path {file_path_input} "
            f"--result_path {output_folder_path}"
        ),
        shell=True,
        check=True,
    )
    subprocess.run(
        (
            "python /TFvelo/TFvelo_analysis_demo.py "
            f"--dataset_path {file_path_input} "
            f"--result_path {output_folder_path}"
        ),
        shell=True,
        check=True,
    )


# Iterates over all sample directories and processes them conserving the directory structure.
for root, dirs, files in os.walk(INPUT_FOLDER):
    for file in files:
        if file.casefold().endswith("clustered.h5ad"):
            file_path_input = os.path.join(root, file)
            output_folder_path = os.path.join(
                MOUNT_PATHS["output"], root.replace(INPUT_FOLDER, "", 1)
            )
            os.makedirs(output_folder_path, exist_ok=True)
            process_data(file_path_input, output_folder_path)
