"""
Create effective genome size estimates using unique-kmers.py.

"""

import subprocess
import re
import json
import os
import sys
import traceback

ref_genome = snakemake.input["ref_genome"]
egs_json = snakemake.output["egs_json"]
log_file = snakemake.log["std_err"]

k_mers = [50, 60, 65, 75, 100, 150, 200]
results = {}

os.makedirs(os.path.dirname(egs_json), exist_ok=True)
os.makedirs(os.path.dirname(log_file), exist_ok=True)

with open(log_file, "a") as logf:
    logf.write("=== Starting effective genome size estimation ===\n")
    logf.write(f"Reference genome: {ref_genome}\n")

    for k in k_mers:
        logf.write(f"--- Running unique-kmers.py with k={k} ---\n")

        cmd = ["unique-kmers.py", "-k", str(k), ref_genome]

        try:
            process = subprocess.run(cmd, capture_output=True, text=True, check=True)

            match = re.search(
                r"Total estimated number of unique \d+-mers:\s+(\d+)", process.stderr
            )
            if match:
                results[str(k)] = int(match.group(1))
            else:
                logf.write("ERROR: Could not parse unique k-mer count\n")
                logf.write("---- STDERR ----\n")
                logf.write(process.stderr + "\n")
                raise RuntimeError(f"Parsing failed for k={k}")

            logf.write(f"SUCCESS: k={k}, unique_kmers={results[str(k)]}\n\n")

        except subprocess.CalledProcessError as e:
            logf.write("ERROR: unique-kmers.py failed\n")
            logf.write("---- STDERR ----\n")
            logf.write((e.stderr or "") + "\n")
            logf.write("---- STDOUT ----\n")
            logf.write((e.stdout or "") + "\n")
            raise

        except Exception:
            logf.write("ERROR: Unexpected Python exception\n")
            logf.write(traceback.format_exc() + "\n")
            raise

    logf.write("=== All k-mer calculations completed successfully ===\n")

with open(egs_json, "w") as f:
    json.dump(results, f, indent=4)
