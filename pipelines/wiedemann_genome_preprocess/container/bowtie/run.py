#!/usr/bin/python
"""This module runs the bowtie2 indexing."""

import json
import math
import multiprocessing
import os
import sys

BASE_COMMAND = "/bowtie2/bowtie2-build"
MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))

# If a specific environment variable is set, appends the respective option.
options = ""

threads = math.floor(multiprocessing.cpu_count() * 0.8)
if threads > 0:
    options += f" --threads {threads}"

if not options:
    print("Running with default options.")
else:
    print("Specified options:" + options)

full_command = f"{BASE_COMMAND}{options} {MOUNT_PATHS['dependencies']['unpack']}/genome.fa {MOUNT_PATHS['output']}/genome"
exit_code = os.waitstatus_to_exitcode(os.system(full_command))
if exit_code != 0:
    sys.exit(exit_code)
