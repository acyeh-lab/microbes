#!/usr/bin/env python3

# Summary:
#   Identify "shared" regions in a chosen reference sequence that are conserved
#   across a directory of FASTA files (each file treated as one strain).
#
#   The script uses exact k-mer matching (default k=35). It:
#     1) Loads all *.fasta files in the given directory (each file = one strain).
#     2) For each strain, builds a UNION set of all k-mers observed anywhere in
#        that file's records (contigs/chromosomes).
#     3) Selects a reference record (sequence) from which all candidate k-mers are taken.
#     4) Counts, for each reference k-mer, in how many strains (files) that k-mer appears.
#     5) Keeps only k-mers that meet a "strain support" threshold:
#          - If --min-fraction is set: require ceil(min_fraction * #strains)
#          - Else if --require-at-least is set: require that many strains
#          - Else: require presence in ALL strains
#     6) Marks base positions on the reference that lie within any supported k-mer,
#        then extracts continuous regions on the reference where:
#          - Region length >= k
#          - Every k-mer window inside the region is supported ("shared")
#     7) Merges overlapping/contiguous regions and writes them to a TSV.
#
# What "shared" means here:
#   A k-mer is "shared" if it is present (exact string match) in at least the
#   required number of STRAINS (files). Presence is counted by strain, not by
#   number of contigs/records.
#
# Inputs (CLI arguments):
#   1) fasta_dir (positional): Directory containing *.fasta files (each file = one strain)
#   2) -k / --kmer: k-mer length, and minimum region length (default 35)
#   3) --min-fraction: fraction of strains required to contain the k-mer (0..1]
#   4) --require-at-least: integer # of strains required (ignored if --min-fraction is set)
#   5) --ref-file: choose a specific reference strain FASTA file (uses its first record)
#   6) --ref-name: substring match to choose a specific reference record name (if --ref-file not set)
#   7) -o / --out: output TSV path (default shared_regions.tsv)
#
# Reference selection (priority order):
#   1) --ref-file: first FASTA record in that file
#   2) --ref-name: first record whose constructed name contains the substring
#   3) default: first record in the first *.fasta file (alphabetical order)
#
# Output:
#   TSV file with a small header block plus tab-separated columns:
#     sequence   start   end
#   Where:
#     - 'sequence' is the reference substring for each shared region
#     - 'start' is 1-based on the reference
#     - 'end' is written as the reference slice end index (see code); interpret as
#       1-based inclusive end if using start=s+1 and end=e for the slice [s:e)
#
# Example:
#   python shared_regions.py genomes_dir/ -k 35 --min-fraction 0.9 --ref-file GCF_000001.fasta -o shared.tsv
#
# Notes / caveats:
#   - Only files ending with *.fasta are read (not .fa/.fna/.gz unless modified).
#   - Reverse complements are NOT considered; matches are strand-sensitive.
#   - "Strain presence" is based on k-mer existence anywhere in that strain file.
#   - Runtime can be high for large reference sequences because support counting
#     scales roughly with (#reference_kmers × #strains).
#   - Memory can be high because per-strain k-mer UNION sets are stored in RAM.


import argparse
import math
from pathlib import Path
from typing import Optional

def read_fasta(path: Path):
    """Yield (header, sequence) for all records in a FASTA file."""
    header = None
    chunks = []
    with path.open("r") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks).upper().replace("U", "T")
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line)
        if header is not None:
            yield header, "".join(chunks).upper().replace("U", "T")


def kmers_positions(seq: str, k: int):
    """Return dict: kmer -> list of start positions (0-based) in seq. Skips kmers containing 'N'."""
    d = {}
    n = len(seq)
    if n < k:
        return d
    for i in range(n - k + 1):
        kmer = seq[i:i+k]
        if 'N' in kmer:
            continue
        d.setdefault(kmer, []).append(i)
    return d


def kmers_set(seq: str, k: int):
    """Return set of k-mers (no positions), skipping kmers containing 'N'."""
    out = set()
    n = len(seq)
    if n < k:
        return out
    for i in range(n - k + 1):
        kmer = seq[i:i+k]
        if 'N' in kmer:
            continue
        out.add(kmer)
    return out


def load_records_and_strain_kmer_sets(fasta_dir: Path, k: int):
    """
    Returns:
      records: list of (strain_id, record_name, sequence) across all .fasta files
      strain_kmers: dict[strain_id] -> set of kmers (union across that file's records)
    """
    records = []
    strain_kmers = {}
    for p in sorted(fasta_dir.glob("*.fasta")):
        strain_id = p.stem  # each file is one strain
        union = set()
        for hdr, seq in read_fasta(p):
            record_name = f"{strain_id}|{hdr.split()[0]}"
            records.append((strain_id, record_name, seq))
            union |= kmers_set(seq, k)
        strain_kmers[strain_id] = union
    return records, strain_kmers


def compute_required(n_strains: int, min_fraction: Optional[float], require_at_least: int):
    """
    Decide how many STRAINS must have the k-mer.
    Priority: --min-fraction (if provided) > --require-at-least (>0) > default (= all strains).
    """
    if min_fraction is not None:
        f = float(min_fraction)
        if f < 0 or f > 1.0:
            raise SystemExit(f"--min-fraction must be in [0,1], got {f}")
        if f == 0.0:
            return 1
        return max(1, math.ceil(f * n_strains))
    if require_at_least and require_at_least > 0:
        return require_at_least
    return n_strains


def main():
    ap = argparse.ArgumentParser(
        description="Find shared continuous sequences (exact matches) of length >= k across FASTA files, "
                    "counting support by number of strains (files), not total records."
    )
    ap.add_argument("fasta_dir", help="Directory containing .fasta files (each file = one strain)")
    ap.add_argument("-k", "--kmer", type=int, default=35, help="Minimum k (default: 35)")
    ap.add_argument("--require-at-least", type=int, default=0,
                    help="Require presence in at least this many strains (default: 0 -> all strains). "
                         "Ignored if --min-fraction is set.")
    ap.add_argument("--min-fraction", type=float, default=None,
                    help="Fraction of strains that must contain a k-mer (0..1]. "
                         "0 => keep k-mers seen in >=1 strain. Takes precedence over --require-at-least.")
    ap.add_argument("--ref-name", type=str, default=None,
                    help="Substring to select a specific reference record for coordinates "
                         "(searches across all records' names). Defaults to the first record.")
    ap.add_argument("--ref-file", type=str, default=None,
                    help="Filename of the reference strain (e.g. GCF_029910175.1.fasta). "
                         "Overrides --ref-name if provided.")
    ap.add_argument("-o", "--out", type=str, default="shared_regions.tsv",
                    help="Output TSV with columns: sequence, start, end (1-based on reference)")
    args = ap.parse_args()

    fasta_dir = Path(args.fasta_dir)
    if not fasta_dir.is_dir():
        raise SystemExit(f"Not a directory: {fasta_dir}")

    # Preload records and per-strain k-mer unions
    k = args.kmer
    records, strain_kmers = load_records_and_strain_kmer_sets(fasta_dir, k)
    if not records:
        raise SystemExit("No FASTA records found.")

    # Pick reference record
    if args.ref_file:
        ref_path = fasta_dir / args.ref_file
        if not ref_path.exists():
            raise SystemExit(f"--ref-file '{args.ref_file}' not found in directory.")
        ref_records = list(read_fasta(ref_path))
        if not ref_records:
            raise SystemExit(f"No FASTA records found in {args.ref_file}.")
        ref_strain = ref_path.stem
        ref_name, ref_seq = ref_records[0]
    elif args.ref_name:
        ref_candidates = [(sid, rname, seq) for (sid, rname, seq) in records if args.ref_name in rname]
        if not ref_candidates:
            raise SystemExit(f"--ref-name '{args.ref_name}' not found among record names.")
        ref_strain, ref_name, ref_seq = ref_candidates[0]
    else:
        ref_strain, ref_name, ref_seq = records[0]

    # Build reference k-mer index
    ref_index = kmers_positions(ref_seq, k)
    if not ref_index:
        raise SystemExit(f"Reference record '{ref_name}' shorter than k={k} or no valid k-mers.")

    # Determine strain-based support threshold
    n_strains = len(strain_kmers)
    required = compute_required(n_strains, args.min_fraction, args.require_at_least)

    # Count in how many STRAINS each reference k-mer appears
    kmer_support = {}
    ref_union = strain_kmers.get(ref_strain, set())
    for km in ref_index.keys():
        kmer_support[km] = 1 if km in ref_union else 0

    for sid, union in strain_kmers.items():
        if sid == ref_strain:
            continue
        for km in kmer_support.keys():
            if km in union:
                kmer_support[km] += 1

    # Keep only reference kmers supported by enough STRAINS
    shared_kmers = {km for km, c in kmer_support.items() if c >= required}
    if not shared_kmers:
        need = (f"ceil({args.min_fraction}*{n_strains})={required}" if args.min_fraction is not None
                else (f"{required}" if args.require_at_least > 0 else "all strains"))
        raise SystemExit(f"No shared k-mers (k={k}) reached threshold across {need}.")

    # Mark shared positions on the reference record
    n = len(ref_seq)
    is_shared = [False] * n
    for km in shared_kmers:
        for pos in ref_index[km]:
            for i in range(pos, pos + k):
                is_shared[i] = True

    # Extract continuous regions
    regions = []
    i = 0
    while i < n:
        while i < n and not is_shared[i]:
            i += 1
        if i >= n:
            break
        j = i
        while j < n and is_shared[j]:
            j += 1
        if j - i >= k:
            ok = True
            for s in range(i, j - k + 1):
                if ref_seq[s:s+k] not in shared_kmers:
                    ok = False
                    break
            if ok:
                regions.append((i, j))
        i = j

    # Merge contiguous/overlapping intervals
    merged = []
    for start, end in regions:
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)

    # Summary stats
    total_bp = sum(end - start for start, end in merged)
    total_records = len(records)
    total_strains = len(strain_kmers)
    used_fraction = args.min_fraction if args.min_fraction is not None else "None"

    # Write output with header block
    with open(args.out, "w") as out:
        out.write(f">Reference_record: {ref_name}\n")
        out.write(f"# Reference_file: {args.ref_file if args.ref_file else 'Auto-selected'}\n")
        out.write(f"# Strains_analyzed: {total_strains}\n")
        out.write(f"# Total_unique_sequences: {total_records}\n")
        out.write(f"# min_fraction_used: {used_fraction}\n")
        out.write(f"# k={k}, Required_support={required}, Total_shared_regions={len(merged)}, Total_bp_matched={total_bp}\n")
        out.write("sequence\tstart\tend\n")
        for s, e in merged:
            subseq = ref_seq[s:e]
            out.write(f"{subseq}\t{s+1}\t{e}\n")

    # Console summary
    print(f"Reference record: {ref_name} (strain={ref_strain})")
    if args.ref_file:
        print(f"Reference file: {args.ref_file}")
    print(f"Strains (files): {total_strains}")
    print(f"Total unique sequences analyzed: {total_records}")
    print(f"min_fraction used: {used_fraction}")
    if args.min_fraction is not None:
        print(f"k = {k}, strain support >= ceil({args.min_fraction} * {total_strains}) = {required}")
    else:
        th = required if args.require_at_least > 0 else total_strains
        print(f"k = {k}, strain support >= {th}")
    print(f"Shared regions: {len(merged)}")
    print(f"Total basepairs matched: {total_bp}")
    print(f"Wrote: {args.out}")


if __name__ == "__main__":
    main()

