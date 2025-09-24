#!/usr/bin/env python3

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
