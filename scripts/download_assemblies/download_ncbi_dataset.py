#!/usr/bin/env python3

import argparse
import subprocess
import os
import zipfile
import sys


def download_dataset(organism, outdir, reference=False):
    import shutil

    os.makedirs(outdir, exist_ok=True)

    print(f"Downloading: {organism}")
    zip_path = os.path.join(outdir, "ncbi_dataset.zip")

    # Build the datasets command
    command = [
        "datasets", "download", "genome", "taxon", organism,
        "--annotated",
        "--assembly-level", "complete",
        "--include", "genome,protein,gff3",
        "--filename", zip_path
    ]

    if reference:
        command.append("--reference")

    # Run download
    try:
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as e:
        if reference:
            print(f"[INFO] No reference genome found for: {organism}. Try running without '--reference'.")
        else:
            print(f"Failed to download data for: {organism}")
        sys.exit(1)

    # Unzip
    print(f"Unzipping: {zip_path}")
    try:
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            zip_ref.extractall(outdir)
        print("Unzip complete")
    except zipfile.BadZipFile:
        print(f"Failed to unzip {zip_path} — file may be corrupt.")
        sys.exit(1)

    # Generate data_summary.tsv
    print("Generating data_summary.tsv...")
    summary_path = os.path.join(outdir, "ncbi_dataset", "data", "data_summary.tsv")

    try:
        result = subprocess.run([
            "dataformat", "tsv", "genome",
            "--package", zip_path,
            "--fields", "organism-name,organism-common-name,organism-infraspecific-strain,organism-tax-id,assminfo-name,accession,source_database,annotinfo-name,assminfo-level,assmstats-contig-n50,assmstats-total-sequence-len,assminfo-biosample-submission-date,annotinfo-featcount-gene-total,assminfo-bioproject,assminfo-biosample-accession"
        ], capture_output=True, check=True, text=True)

        with open(summary_path, "w") as f:
            f.write(result.stdout)
        print(f"Saved data_summary.tsv to {summary_path}")

    except subprocess.CalledProcessError as e:
        print("Failed to generate data_summary.tsv")
        print(e.stderr)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download NCBI genome data for a given organism.")
    parser.add_argument("--organism", required=True, help="Scientific name of the organism (in quotes)")
    parser.add_argument("--outdir", required=True, help="Directory to save the downloaded and unzipped data")
    parser.add_argument("--reference", action="store_true", help="Only download reference genomes (optional)")
    args = parser.parse_args()

    download_dataset(args.organism, args.outdir, args.reference)

