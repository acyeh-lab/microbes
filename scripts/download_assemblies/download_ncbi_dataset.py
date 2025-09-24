import argparse
import subprocess
import os
import zipfile
import sys

def download_dataset(organism, outdir):
    os.makedirs(outdir, exist_ok=True)

    print(f"Downloading: {organism}")
    zip_path = os.path.join(outdir, "ncbi_dataset.zip")

    try:
        subprocess.run([
            "datasets", "download", "genome", "taxon", organism,
            "--annotated",
            "--assembly-level", "complete",
            "--include", "genome,protein,gff3",
            "--filename", zip_path
        ], check=True)
    except subprocess.CalledProcessError:
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

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download NCBI genome data for a given organism.")
    parser.add_argument("--organism", required=True, help="Scientific name of the organism (in quotes)")
    parser.add_argument("--outdir", required=True, help="Directory to save the downloaded and unzipped data")
    args = parser.parse_args()

    download_dataset(args.organism, args.outdir)

