#!/usr/bin/python
"""This module downloads blacklists."""

import json
import logging
import os
import shutil
import sys

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
TEMPORARY_OUTPUT_PATH = os.path.join(MOUNT_PATHS["output"], "tmp")

# Setup.
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
os.makedirs(TEMPORARY_OUTPUT_PATH, exist_ok=True)

logging.info("Downloading blacklist...")
command_download_blacklist = (
    f"git -C {TEMPORARY_OUTPUT_PATH} clone https://github.com/Boyle-Lab/Blacklist.git && "
    f"git -C {os.path.join(TEMPORARY_OUTPUT_PATH, 'Blacklist')} checkout 61a04d2"
)
exit_code_download_blacklist = os.waitstatus_to_exitcode(
    os.system(command_download_blacklist)
)
if exit_code_download_blacklist != 0:
    sys.exit(exit_code_download_blacklist)

logging.info("Extracting blacklist...")
command_unzip_blacklist = f"gzip -d {os.path.join(TEMPORARY_OUTPUT_PATH, 'Blacklist', 'lists', 'mm10-blacklist.v2.bed.gz')}"
exit_code_unzip_blacklist = os.waitstatus_to_exitcode(
    os.system(command_unzip_blacklist)
)
if exit_code_unzip_blacklist != 0:
    sys.exit(exit_code_unzip_blacklist)

logging.info("Renaming files...")
os.rename(
    os.path.join(TEMPORARY_OUTPUT_PATH, "Blacklist", "lists", "mm10-blacklist.v2.bed"),
    os.path.join(MOUNT_PATHS["output"], "blacklist.bed"),
)

logging.info("Cleaning up temporary files...")
shutil.rmtree(TEMPORARY_OUTPUT_PATH)
logging.info("Done.")
