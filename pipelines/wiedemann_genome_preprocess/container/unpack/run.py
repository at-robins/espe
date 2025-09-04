#!/usr/bin/python
"""This module downloads the genome files."""

import json
import logging
import os
import sys

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
DOWNLOAD_PATH_GENOME = os.path.join(MOUNT_PATHS["output"], "genome.fa.gz")
DOWNLOAD_PATH_ANNOTATION_GTF = os.path.join(MOUNT_PATHS["output"], "annotations.gtf.gz")
DOWNLOAD_PATH_ANNOTATION_GFF = os.path.join(MOUNT_PATHS["output"], "annotations.gff.gz")

# Setup.
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
os.makedirs(MOUNT_PATHS["output"], exist_ok=True)


def annotation_url(organism, release, gtf):
    """
    Returns the URL to download the specific annotation file.
    """
    if gtf:
        extension_name = "gtf"
    else:
        extension_name = "gff3"

    return (
        f"https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_{organism}"
        f"/release_{release}/gencode.v{release}.annotation.{extension_name}.gz"
    )


organism_env = os.environ.get("GLOBAL_ORGANISM")
if organism_env is not None and organism_env == "human":
    print("\tUsing human data (GRCh38.p14)...", flush=True)
    genome_url = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/p14/hg38.p14.fa.gz"
    )
    annotation_release_id = "48"
    annotation_organism_id = "human"
else:
    print("\tUsing mouse data (GRCm38.p6)...", flush=True)
    genome_url = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/p6/mm10.p6.fa.gz"
    )
    annotation_release_id = "M25"
    annotation_organism_id = "mouse"

logging.info("Downloading genome...")
command_download_genome = (
    "wget " "--no-verbose " f"-O {DOWNLOAD_PATH_GENOME} " f"{genome_url}"
)
exit_code_download_genome = os.waitstatus_to_exitcode(
    os.system(command_download_genome)
)
if exit_code_download_genome != 0:
    sys.exit(exit_code_download_genome)

logging.info("Extracting genome...")
command_unzip_genome = f"gunzip {DOWNLOAD_PATH_GENOME}"
exit_code_unzip_genome = os.waitstatus_to_exitcode(os.system(command_unzip_genome))
if exit_code_unzip_genome != 0:
    sys.exit(exit_code_unzip_genome)


logging.info("Downloading annotations...")
command_download_annotations = (
    "wget "
    "--no-verbose "
    f"-O {DOWNLOAD_PATH_ANNOTATION_GTF} "
    f"{annotation_url(organism=annotation_organism_id, release=annotation_release_id, gtf=True)} && "
    "wget "
    "--no-verbose "
    f"-O {DOWNLOAD_PATH_ANNOTATION_GFF} "
    f"{annotation_url(organism=annotation_organism_id, release=annotation_release_id, gtf=False)}"
)
exit_code_download_annotations = os.waitstatus_to_exitcode(
    os.system(command_download_annotations)
)
if exit_code_download_annotations != 0:
    sys.exit(exit_code_download_annotations)

logging.info("Extracting annotations...")
command_unzip_annotations = (
    f"gunzip {DOWNLOAD_PATH_ANNOTATION_GTF} && "
    f"gunzip {DOWNLOAD_PATH_ANNOTATION_GFF}"
)
exit_code_unzip_annotations = os.waitstatus_to_exitcode(
    os.system(command_unzip_annotations)
)
if exit_code_unzip_annotations != 0:
    sys.exit(exit_code_unzip_annotations)

logging.info("Done.")
