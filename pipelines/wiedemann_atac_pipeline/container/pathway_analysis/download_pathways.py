#!/usr/bin/python
"""This module downloads the pathway networks."""

print("Loading decoupler...", flush=True)
import decoupler
import os

OUTPUT_PATH = "/pathway_networks"
os.makedirs(OUTPUT_PATH, exist_ok=True)
print("Downloading human pathway networks...", flush=True)
decoupler.get_progeny(organism="human").to_pickle(
    os.path.join(OUTPUT_PATH, "human.pkl")
)
print("Downloading murine pathway networks...", flush=True)
decoupler.get_progeny(organism="mouse").to_pickle(
    os.path.join(OUTPUT_PATH, "mouse.pkl")
)
print("Done!", flush=True)