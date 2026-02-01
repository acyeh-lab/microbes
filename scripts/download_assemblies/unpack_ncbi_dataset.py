#!/usr/bin/env python3

# Summary:
#   Unpack genome and protein FASTA files from an NCBI Datasets genome package
#   (the "ncbi_dataset/" directory) and generate per-folder sample sheets that
#   combine file paths with NCBI assembly metadata.
#
#   The script:
#     - Reads an NCBI Datasets "manifest" (a JSON structure expected to contain
#       a top-level key 'assemblies').
#     - For each assembly in manifest['assemblies']:
#         * Locates files of type GENOMIC_NUCLEOTIDE_FASTA and PROTEIN_FASTA
#         * Copies each file from ncbi_dataset/data/<filePath> into:
#             - genomes/   (for genomic FASTA)
#             - proteins/  (for protein FASTA)
#         * Writes the copied file as gzip-compressed output (*.gz), even if the
#           source file is not gzipped.
#         * Builds two in-memory manifest lists:
#             - genomic_manifest: [{sample: <accession>, file: <filename>}, ...]
#             - protein_manifest: [{sample: <accession>, file: <filename>}, ...]
#     - Writes two CSV sample sheets:
#         - genomes/samplesheet.csv
#         - proteins/samplesheet.csv
#       Each samplesheet includes:
#         - sample (assembly accession)
#         - file (copied filename in genomes/ or proteins/)
#       And merges in metadata from:
#         - ncbi_dataset/data/data_summary.tsv
#       via Assembly Accession.
#
# Inputs / expected files:
#   - An NCBI Datasets output directory structure in the current working directory:
#       ncbi_dataset/
#         data/
#           data_summary.tsv
#           <files referenced by manifest.json>
#   - A JSON manifest file (read by read_manifest()) containing:
#       {"assemblies": [{"accession": "...", "files": [...]} , ...]}
#     NOTE: read_manifest() is referenced but not shown in this snippet.
#
# Output:
#   - genomes/
#       *.fna.gz (or similar genome FASTA filenames)
#       samplesheet.csv
#   - proteins/
#       *_protein.faa.gz
#       samplesheet.csv
#
# Example usage:
#   python unpack_ncbi_dataset_to_samplesheets.py
#   (Run in a directory that contains ncbi_dataset/ and the manifest JSON.)
#
# Notes / caveats:
#   - This script assumes the NCBI metadata table has a column named
#     "Assembly Accession" and merges on that field.
#   - The genome filename is derived from the original filePath basename.
#     The protein filename is forced to "<accession>_protein.faa.gz".
#   - copy_file() reads the entire input file into memory (txt = open(...).read()).
#     For large FASTA files this can be very memory-inefficient; streaming copy
#     is preferable.
#   - The provided snippet ends mid-line ("hand..."), suggesting copy_file() is
#     incomplete as pasted; as written it would raise an exception.
#   - The code appends ".gz" to outputs unconditionally; it does not preserve
#     original compression state.



import gzip
import json
import pandas as pd
from pathlib import Path

def main():
    manifest = read_manifest()

    assert "assemblies" in manifest, "Expected to find 'assemblies' in manifest"

    genomic_manifest = []
    protein_manifest = []

    for assembly in manifest['assemblies']:
        unpack_assembly(assembly, genomic_manifest, protein_manifest)

    write_manifest(genomic_manifest, "genomes/samplesheet.csv")
    write_manifest(protein_manifest, "proteins/samplesheet.csv")


def write_manifest(manifest, fpo):
    manifest = pd.DataFrame(manifest)

    # Merge in the metadata
    meta = pd.read_csv("ncbi_dataset/data/data_summary.tsv", sep="\t")

    manifest = manifest.merge(meta, left_on="sample", right_on="Assembly Accession")
    manifest.to_csv(fpo, index=None)


def unpack_assembly(assembly: dict, genomic_manifest: list, protein_manifest: list):
    for kw in ['accession', 'files']:
        if kw not in assembly:
            return

    print(f"Processing {assembly['accession']}")

    for file in assembly['files']:
        if file['fileType'] == 'GENOMIC_NUCLEOTIDE_FASTA':
            fn = file['filePath'].split("/")[-1]
            if not fn.endswith(".gz"):
                fn = fn + ".gz"
            copy_file("ncbi_dataset/data/" + file['filePath'], f"genomes/{fn}")
            genomic_manifest.append(dict(sample=assembly['accession'], file=fn))

        elif file['fileType'] == 'PROTEIN_FASTA':
            fn = f"{assembly['accession']}_protein.faa"
            if not fn.endswith(".gz"):
                fn = fn + ".gz"
            copy_file("ncbi_dataset/data/" + file['filePath'], f"proteins/{fn}")
            protein_manifest.append(dict(sample=assembly['accession'], file=fn))


def copy_file(file_in: str, file_out: str):
    txt = open(file_in).read()
    if not file_out.endswith('.gz'):
        file_out = file_out + '.gz'
    file_out = Path(file_out)
    file_out.parent.mkdir(exist_ok=True, parents=True)
    with gzip.open(file_out, 'wt') as handle:
        handle.write(txt)


def read_manifest(fp="ncbi_dataset/data/dataset_catalog.json"):

    file = Path(fp)
    assert file.exists(), f"Expected file to exist: {fp}"

    with open(file) as handle:
        return json.load(handle)


if __name__ == "__main__":
    main()
