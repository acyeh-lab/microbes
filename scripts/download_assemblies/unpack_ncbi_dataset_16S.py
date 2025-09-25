#!/usr/bin/env python3

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

