#!/usr/bin/python
"""This module determines the mappable genome size."""

import json
import os
import shutil

from bs4 import BeautifulSoup

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER = next(iter(MOUNT_PATHS["dependencies"].values()))
GENOME_PATH = os.path.join(MOUNT_PATHS["globals"]["GENOME"], "genome.fa")
SEQUENCE_LENGTH_TAG = "<td>Sequence length</td>"

kmer_map = {}

# Iterates over all sample directories and processes them conserving the directory structure.
for root, dirs, files in os.walk(INPUT_FOLDER):
    read_lengths = []
    for file in files:
        if file.casefold().endswith(".html"):
            file_input_path = os.path.join(root, file)
            print(f"Processing file {file_input_path}...", flush=True)

            # Parses the FastQC HTML report.
            with open(file=file_input_path, mode="rt", encoding="utf-8") as qc_file:
                qc_html = BeautifulSoup(qc_file, "html.parser")
                sequence_length_string = qc_html.find(
                    name="td", string="Sequence length"
                ).next_sibling.string
                # Read length is specified as range in the format "[min]-[max]"
                read_lengths.extend(
                    [int(x.strip()) for x in sequence_length_string.split("-")]
                )

    if len(read_lengths) > 0:
        folder_output_path = os.path.join(
            MOUNT_PATHS["output"],
            os.path.normpath(os.path.relpath(root, INPUT_FOLDER)),
        )
        os.makedirs(folder_output_path, exist_ok=True)
        report_file_path = os.path.join(folder_output_path, "unique_kmer_report.txt")
        max_read_length = max(read_lengths)
        cached_report = kmer_map.get(max_read_length)
        
        if cached_report is not None:
            print(f"\tFound cached kmers for size {max_read_length}...", flush=True)
            shutil.copy(cached_report, report_file_path)
        else:
            print(f"\tCalculating kmers for size {max_read_length}...", flush=True)
            
            full_command = (
                f"unique-kmers.py "
                f"--ksize {max_read_length} "
                f"--report {report_file_path} "
                f"{GENOME_PATH}"
            )
            exit_code = os.waitstatus_to_exitcode(os.system(full_command))
            if exit_code != 0:
                sys.exit(exit_code)

            # Caches the result as this is will be the same for equal read lengths.
            kmer_map[max_read_length] = report_file_path