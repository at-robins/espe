#!/usr/bin/python
"""This module runs the bowtie2 indexing."""

import json
import logging
import math
import multiprocessing
import os
import sys

BASE_COMMAND = "bowtie2-build"
MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))

# Setup.
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")

# If a specific environment variable is set, appends the respective option.
options = ""

threads = math.floor(multiprocessing.cpu_count() * 0.8)
if threads > 0:
    options += f" --threads {threads}"

if not options:
    logging.info("Running with default options.")
else:
    logging.info(f"Specified options: {options}")

logging.info("Creating indices...")
full_command = f"{BASE_COMMAND}{options} {MOUNT_PATHS['dependencies']['get_genome']}/genome.fa {MOUNT_PATHS['output']}/genome"
exit_code = os.waitstatus_to_exitcode(os.system(full_command))
if exit_code != 0:
    sys.exit(exit_code)

logging.info("Done.")
