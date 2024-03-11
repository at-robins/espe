#!/usr/bin/python
"""This module downloads the genome files."""

import json
import logging
import os
import shutil
import sys

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
TEMPORARY_OUTPUT_PATH = os.path.join(MOUNT_PATHS["output"], "tmp")
DOWNLOAD_PATH_GENOME = os.path.join(TEMPORARY_OUTPUT_PATH, "genome.zip")

# Setup.
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
os.makedirs(TEMPORARY_OUTPUT_PATH, exist_ok=True)

logging.info("Downloading genome and annotations...")
command_download_genome = (
    "wget "
    "--no-verbose "
    f"-O {DOWNLOAD_PATH_GENOME} "
    "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_000001635.26/download"
    "?include_annotation_type=GENOME_FASTA,GENOME_GFF"
)
exit_code_download_genome = os.waitstatus_to_exitcode(
    os.system(command_download_genome)
)
if exit_code_download_genome != 0:
    sys.exit(exit_code_download_genome)

logging.info("Extracting genome and annotations...")
command_unzip_genome = f"unzip -d {TEMPORARY_OUTPUT_PATH} {DOWNLOAD_PATH_GENOME}"
exit_code_unzip_genome = os.waitstatus_to_exitcode(os.system(command_unzip_genome))
if exit_code_unzip_genome != 0:
    sys.exit(exit_code_unzip_genome)

logging.info("Renaming files...")
os.rename(
    os.path.join(
        TEMPORARY_OUTPUT_PATH,
        "ncbi_dataset",
        "data",
        "GCF_000001635.26",
        "GCF_000001635.26_GRCm38.p6_genomic.fna",
    ),
    os.path.join(MOUNT_PATHS["output"], "genome.fa"),
)
os.rename(
    os.path.join(
        TEMPORARY_OUTPUT_PATH,
        "ncbi_dataset",
        "data",
        "GCF_000001635.26",
        "genomic.gff",
    ),
    os.path.join(MOUNT_PATHS["output"], "annotations.gff"),
)

logging.info("Compressing files...")
command_zip = f"gzip -k {os.path.join(MOUNT_PATHS['output'], 'genome.fa')}"
exit_code_zip = os.waitstatus_to_exitcode(os.system(command_zip))
if exit_code_zip != 0:
    sys.exit(exit_code_zip)

logging.info("Cleaning up temporary files...")
shutil.rmtree(TEMPORARY_OUTPUT_PATH)
logging.info("Done.")
