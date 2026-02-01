#!/usr/bin/env python3
# Summary:
#   Download and unpack NCBI genome datasets for a specified organism using the
#   NCBI Datasets CLI, with optional restriction to reference genomes only.
#
#   The script:
#     - Invokes the `datasets` command-line tool to download complete, annotated
#       genome assemblies for a given taxon (scientific name).
#     - Requests genome FASTA, protein FASTA, and GFF3 annotation files.
#     - Restricts downloads to assemblies at the "complete" level.
#     - Optionally limits results to reference genomes via the --reference flag.
#     - Verifies that the downloaded file is a valid ZIP archive.
#     - Unzips the dataset into the specified output directory.
#     - Generates a tab-delimited metadata file (data_summary.tsv) using the
#       `dataformat` CLI, capturing key assembly and annotation statistics.
#
# Inputs (CLI arguments):
#   --organism   Scientific name of the organism or taxon (e.g. "Escherichia coli")
#   --outdir     Output directory where data will be downloaded and extracted
#   --reference  (Optional) Restrict downloads to reference genomes only
#
# External dependencies:
#   - NCBI Datasets CLI (`datasets`)
#   - NCBI Dataformat CLI (`dataformat`)
#   Both tools must be installed and available on $PATH.
#
# Output:
#   Within --outdir:
#     - ncbi_dataset.zip
#         ZIP archive produced by the Datasets CLI
#     - ncbi_dataset/
#         Extracted directory tree containing genome FASTA, protein FASTA,
#         GFF3 annotation files, and metadata
#     - ncbi_dataset/data/data_summary.tsv
#         Tab-separated summary of assemblies, including:
#           organism name, strain, tax_id, accession, assembly level,
#           N50, genome length, gene count, BioProject, and BioSample
#
# Example:
#   python download_ncbi_genomes.py \
#       --organism "Escherichia coli" \
#       --outdir ./ecoli_genomes
#
#   python download_ncbi_genomes.py \
#       --organism "Escherichia coli" \
#       --outdir ./ecoli_reference_genomes \
#       --reference
#
# Notes / caveats:
#   - The script skips organisms for which no reference genome exists when
#     --reference is specified (does not exit the program).
#   - Only assemblies with assembly level "complete" are requested; draft or
#     chromosome-level assemblies are excluded.
#   - The script assumes a working internet connection and sufficient disk space.
#   - The ZIP file is removed only if detected as invalid; otherwise it is retained.
#   - Errors during download, unzip, or metadata generation are reported but do
#     not cause a hard exit, allowing batch-style usage.


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
    except subprocess.CalledProcessError:
        if reference:
            print(f"[INFO] No reference genome found for: {organism}. Try running without '--reference'.")
        else:
            print(f"[ERROR] Failed to download data for: {organism}")
        return  # Skip instead of exit

    # Check that the file exists and is a valid zip
    if not os.path.exists(zip_path):
        print(f"[ERROR] ZIP file not created for {organism}. Skipping.")
        return

    if not zipfile.is_zipfile(zip_path):
        print(f"[ERROR] ZIP file is invalid or corrupt for {organism}. Skipping.")
        os.remove(zip_path)  # Optional: clean up bad ZIP
        return

    # Unzip
    print(f"Unzipping: {zip_path}")
    try:
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            zip_ref.extractall(outdir)
        print("Unzip complete")
    except zipfile.BadZipFile:
        print(f"[ERROR] Failed to unzip {zip_path} — file may be corrupt. Skipping.")
        return

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
        print(f"[ERROR] Failed to generate data_summary.tsv for {organism}")
        print(e.stderr)
        return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download NCBI genome data for a given organism.")
    parser.add_argument("--organism", required=True, help="Scientific name of the organism (in quotes)")
    parser.add_argument("--outdir", required=True, help="Directory to save the downloaded and unzipped data")
    parser.add_argument("--reference", action="store_true", help="Only download reference genomes (optional)")
    args = parser.parse_args()

    download_dataset(args.organism, args.outdir, args.reference)

