#!/usr/bin/python
"""This module calculates RNA velocity."""

import json
import math
import multiprocessing
import os
import subprocess

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER = next(iter(MOUNT_PATHS["dependencies"].values()))
TFVELO_BASE_PATH = "/TFvelo"
TFVELO_RUN_PATH = os.path.join(TFVELO_BASE_PATH, "TFvelo_run_demo.py")
TFVELO_ANALYSIS_PATH = os.path.join(TFVELO_BASE_PATH, "TFvelo_analysis_demo.py")

threads = math.floor(multiprocessing.cpu_count() * 0.8)
if threads < 1:
    threads = 1

def process_data(file_path_input, output_folder_path):
    """
    Calculates RNA velocity.
    """
    print(f"Processing file {file_path_input}", flush=True)
    subprocess.run(
        (
            f"python {TFVELO_RUN_PATH} "
            f"--n_jobs {threads} "
            f"--dataset_path {file_path_input} "
            f"--result_path {output_folder_path}"
        ),
        cwd=TFVELO_BASE_PATH,
        shell=True,
        check=True,
    )
    subprocess.run(
        (
            f"python {TFVELO_ANALYSIS_PATH} "
            f"--dataset_path {file_path_input} "
            f"--result_path {output_folder_path}"
        ),
        cwd=TFVELO_BASE_PATH,
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
