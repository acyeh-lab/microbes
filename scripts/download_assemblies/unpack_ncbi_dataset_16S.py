#!/usr/bin/env python3
# Summary:
#   Unpack genomic/protein FASTA files (and GFF3 annotations) from an NCBI Datasets
#   genome package, generate per-folder sample sheets with metadata, and extract
#   annotated 16S/23S rRNA sequences directly from the assembly using GFF3 coordinates.
#
#   The script:
#     1) Reads an NCBI Datasets JSON manifest (expected to contain a top-level
#        key 'assemblies').
#     2) For each assembly:
#         - Copies the GENOMIC_NUCLEOTIDE_FASTA into genomes/ as gzipped output
#           and records it in a genome samplesheet manifest.
#         - Copies the PROTEIN_FASTA into proteins/ as gzipped output and records
#           it in a protein samplesheet manifest.
#         - Copies the GFF3 annotation file into a standardized per-assembly path:
#             ncbi_dataset/data/<accession>/genomic.gff
#           (written uncompressed when compress=False).
#         - If both the GFF and FASTA are available, parses the GFF to locate rRNA
#           features annotated as 16S or 23S and extracts those sequences from the
#           genome FASTA, writing them to:
#             16S_rRNA/<accession>.fasta
#             23S_rRNA/<accession>.fasta
#           Respecting the strand (reverse-complement on '-' strand).
#     3) Writes two merged sample sheets:
#         - genomes/samplesheet.csv
#         - proteins/samplesheet.csv
#        by merging (sample,file) with NCBI metadata from data_summary.tsv
#        on "Assembly Accession".
#
# Inputs / expected files:
#   - NCBI Datasets output directory in the current working directory:
#       ncbi_dataset/
#         data/
#           data_summary.tsv
#           <assembly subfolders and files referenced by the manifest>
#   - A JSON manifest file (read by read_manifest()) with the structure:
#       {"assemblies": [{"accession": "...", "files": [...]} , ...]}
#     NOTE: read_manifest() is referenced but not shown in this snippet.
#
# External dependencies:
#   - pandas
#   - Biopython (Bio.SeqIO)
#
# Outputs:
#   - genomes/
#       <genomic_fasta_basename>.gz
#       samplesheet.csv
#   - proteins/
#       <accession>_protein.faa.gz
#       samplesheet.csv
#   - ncbi_dataset/data/<accession>/genomic.gff     (copied GFF3 per assembly)
#   - 16S_rRNA/
#       <accession>.fasta      (one or more extracted 16S entries)
#   - 23S_rRNA/
#       <accession>.fasta      (one or more extracted 23S entries)
#
# What "extract_rRNA_inline" does:
#   - Loads the genome FASTA into memory as a dict of contig_id → sequence.
#   - Scans the GFF file for rRNA features whose attributes mention "16S" or "23S".
#   - Uses the (seq_id, start, end, strand) fields to slice the corresponding contig.
#   - Reverse-complements the slice for minus-strand features.
#   - Writes each extracted feature as a FASTA record with informative headers.
#
# Example usage:
#   python unpack_and_extract_rrna.py
#   (Run in a directory that contains ncbi_dataset/ and the manifest JSON.)
#
# Notes / caveats:
#   - This script assumes the metadata TSV has a column named "Assembly Accession".
#   - rRNA detection is based on substring checks for "16S" or "23S" in the GFF
#     attributes; annotation conventions vary across genomes, so some rRNAs may be
#     missed if labeled differently.
#   - SeqIO.to_dict loads the entire genome FASTA into RAM; for very large genomes,
#     this can be memory-intensive.
#   - The code standardizes the GFF name to "genomic.gff" and places it under an
#     accession-specific folder; this may overwrite existing files on re-run.
#   - copy_file() is declared but truncated here; as pasted, the script is incomplete
#     until copy_file() (and read_manifest()) are fully implemented.


import gzip
import json
import pandas as pd
from pathlib import Path
from Bio import SeqIO

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
    meta = pd.read_csv("ncbi_dataset/data/data_summary.tsv", sep="\t")
    manifest = manifest.merge(meta, left_on="sample", right_on="Assembly Accession")
    manifest.to_csv(fpo, index=None)


def unpack_assembly(assembly: dict, genomic_manifest: list, protein_manifest: list):
    if "accession" not in assembly or "files" not in assembly:
        return

    accession = assembly["accession"]
    gff_file = None
    fna_file = None

    print(f"Processing {accession}")

    for file in assembly['files']:
        file_path = "ncbi_dataset/data/" + file['filePath']
        file_name = Path(file['filePath']).name

        if file['fileType'] == 'GENOMIC_NUCLEOTIDE_FASTA':
            fn = file_name + ".gz"
            copy_file(file_path, f"genomes/{fn}")
            genomic_manifest.append(dict(sample=accession, file=fn))
            fna_file = Path(file_path)

        elif file['fileType'] == 'PROTEIN_FASTA':
            fn = f"{accession}_protein.faa.gz"
            copy_file(file_path, f"proteins/{fn}")
            protein_manifest.append(dict(sample=accession, file=fn))

        elif file['fileType'] == 'GFF3':
            fn = "genomic.gff"
            dest_path = f"ncbi_dataset/data/{accession}/{fn}"
            copy_file(file_path, dest_path, compress=False)
            gff_file = Path(dest_path)

    if gff_file and fna_file:
        extract_rRNA_inline(gff_file, fna_file, accession)
    else:
        print(f"Missing GFF or FASTA for {accession}, skipping 16S")

def extract_rRNA_inline(gff_path: Path, fna_path: Path, accession: str):
    if not gff_path.exists() or not fna_path.exists():
        print(f"Missing file for {accession}, skipping rRNA")
        return

    seqs = SeqIO.to_dict(SeqIO.parse(open(fna_path), "fasta"))

    output_dirs = {
        "16S": Path("16S_rRNA"),
        "23S": Path("23S_rRNA")
    }
    for path in output_dirs.values():
        path.mkdir(exist_ok=True)

    handles = {
        "16S": open(output_dirs["16S"] / f"{accession}.fasta", "w"),
        "23S": open(output_dirs["23S"] / f"{accession}.fasta", "w")
    }

    found = {"16S": 0, "23S": 0}
    with open(gff_path) as gff_handle:
        for line in gff_handle:
            if '\trRNA\t' in line and ('16S' in line or '23S' in line):
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue

                seq_id, source, feature_type, start, end, score, strand, phase, attributes = parts
                start = int(start) - 1
                end = int(end)
                seq_record = seqs.get(seq_id)

                rna_type = None
                if "16S" in attributes:
                    rna_type = "16S"
                elif "23S" in attributes:
                    rna_type = "23S"

                if rna_type and seq_record:
                    subseq = seq_record.seq[start:end]
                    if strand == '-':
                        subseq = subseq.reverse_complement()

                    found[rna_type] += 1
                    handles[rna_type].write(
                        f">{accession}|{seq_id}:{start+1}-{end}({strand})|{rna_type}_rRNA_{found[rna_type]}\n{subseq}\n"
                    )

    for f in handles.values():
        f.close()

    if sum(found.values()) == 0:
        print(f"No 16S or 23S rRNA found for {accession}")
    else:
        print(f"Found {found['16S']}x 16S and {found['23S']}x 23S rRNA entries for {accession}")

def copy_file(file_in: str, file_out: str, compress: bool = True):
    file_out = Path(file_out)
    file_out.parent.mkdir(exist_ok=True, parents=True)

    with open(file_in, "r") as f:
        txt = f.read()

    if compress:
        if not str(file_out).endswith(".gz"):
            file_out = file_out.with_suffix(file_out.suffix + ".gz")
        with gzip.open(file_out, "wt") as handle:
            handle.write(txt)
    else:
        with open(file_out, "w") as handle:
            handle.write(txt)

def read_manifest(fp="ncbi_dataset/data/dataset_catalog.json"):
    file = Path(fp)
    assert file.exists(), f"Expected file to exist: {fp}"
    with open(file) as handle:
        return json.load(handle)


if __name__ == "__main__":
    main()

