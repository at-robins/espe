#!/usr/bin/python
"""This module calculates the mappable genome size."""

import csv
import json
import logging
import os
import subprocess

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))

# Setup.
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
counts_output_path = os.path.join(MOUNT_PATHS["output"], "genome_counts.tsv")
size_ouput_path = os.path.join(MOUNT_PATHS["output"], "mappable_genome_size.txt")

logging.info("Calculating mappable genome size...")
subprocess.run(
    (
        "faCount "
        f"{MOUNT_PATHS['dependencies']['get_genome']}/genome.fa "
        f"-summary > {counts_output_path}"
    ),
    shell=True,
    check=True,
)

with open(counts_output_path, newline="", encoding="utf-8") as tsvfile:
    info_reader = csv.DictReader(tsvfile, dialect="unix", delimiter="\t", quotechar='"')
    for row in info_reader:
        if row["#seq"] == "total":
            with open(size_ouput_path, mode="wt", encoding="utf-8") as output_file:
                total_size = int(row["A"]) + int(row["T"]) + int(row["G"]) + int(row["C"])
                logging.info(f"Mappable genome size: {total_size} bp")
                output_file.write(str(total_size))
            break

logging.info("Done.")
