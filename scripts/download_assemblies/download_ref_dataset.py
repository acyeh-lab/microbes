#!/usr/bin/env python3

import pandas as pd
import subprocess
import os
import argparse

def download_ncbi_package(taxon_id, sci_name, outdir):
    os.makedirs(outdir, exist_ok=True)
    zip_path = os.path.join(outdir, "ncbi_dataset.zip")

    print(f"Downloading data for {sci_name} (TaxID: {taxon_id})")

    # First attempt: representative/reference genomes
    try:
        subprocess.run([
            "datasets", "download", "genome", "taxon", str(taxon_id),
            "--reference",
            "--annotated",
            ##"--assembly-level", "complete",
            "--include", "genome,protein,gff3",
            "--filename", zip_path
        ], check=True)
    except subprocess.CalledProcessError:
        print(f"No representative genome for {sci_name}. Trying best annotated...")

        # Fallback attempt: best annotated complete genome
        try:
            subprocess.run([
                "datasets", "download", "genome", "taxon", str(taxon_id),
                "--annotated",
                "--assembly-level", "complete",
                "--include", "genome,protein,gff3",
                "--filename", zip_path
            ], check=True)
        except subprocess.CalledProcessError:
            print(f"Failed to download data for {sci_name} (TaxID: {taxon_id})")
            return

    # Unzip result
    try:
        subprocess.run(["unzip", "-o", zip_path, "-d", outdir], check=True)
    except subprocess.CalledProcessError:
        print(f"Failed to unzip for {sci_name} (TaxID: {taxon_id})")

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


def main():
    parser = argparse.ArgumentParser(description="Download NCBI genome data for taxon IDs listed in a CSV.")
    parser.add_argument("--csv", required=True, help="Path to input CSV with columns: NCBI.taxon.ID and scientific.name")
    parser.add_argument("--outdir", default="Downloads", help="Output directory for downloaded data")

    args = parser.parse_args()
    df = pd.read_csv(args.csv)

    for _, row in df.iterrows():
        taxid = row["NCBI.taxon.ID"]
        name = row["scientific.name"].replace(" ", "_")
        outdir = os.path.join(args.outdir, name)
        download_ncbi_package(taxid, name, outdir)

if __name__ == "__main__":
    main()

