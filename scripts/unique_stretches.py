#!/usr/bin/env python3

# Summary:
#   Identify maximal contiguous "unique" stretches in a query sequence (RNA or DNA)
#   relative to a reference database of FASTA files, using exact k-mer presence/absence
#   (default k=35). Additionally, record *where* query k-mers are found in the reference,
#   including the reference FASTA header (the ">" line token) and coordinates.
#
#   This script performs two related tasks:
#
#     (A) Uniqueness scanning (query-centric)
#       - Slide a k-mer window across the query
#       - Mark each query k-mer as present/absent in the reference database
#       - Convert consecutive runs of "absent" k-mers into maximal base-interval stretches
#       - Report stretches with length >= --min_len as candidate "unique" regions
#
#     (B) Hit localization (reference-context reporting)
#       - Build a reference k-mer index that retains provenance:
#           kmer -> list of hits (ref_file, ref_header, ref_pos0, strand)
#       - For each query k-mer that is present, write where it occurs:
#           which reference file, which FASTA record header, start position, and strand
#
# Inputs (CLI arguments):
#   --query        Query sequence as a literal string (RNA or DNA). RNA 'U' is converted to 'T'.
#   --query_fasta  FASTA file containing query sequence(s). Only the FIRST record is used.
#   --ref_dir      Directory containing reference FASTA-like files, searched recursively.
#                 Supported extensions: .fa, .fasta, .fna, .ffn, .fas (also .gz-compressed).
#   -k             k-mer length that defines exact matching (default: 35).
#   --min_len      Minimum length of reported unique stretches (default: 35; forced to >= k).
#   --out_prefix   Output filename prefix (default: unique_stretches). If no path is included,
#                 files are written to the current working directory.
#
# Reference parsing / indexing behavior:
#   - Reference FASTA files can be plain text or gzipped (*.gz). Files are streamed record-by-record
#     (the FASTA parser does not load entire files into memory).
#   - The k-mer index stores BOTH:
#       * the k-mer itself (strand '+')
#       * the reverse complement as a key (strand '-')
#     This makes query matching strand-agnostic without having to also reverse-complement the query.
#   - K-mers containing 'N' or non-ACGT letters are skipped in both reference and query.
#
# Output:
#   1) <out_prefix>.tsv
#        Unique stretches in the query:
#          query_id   start0   end0_excl   length   sequence
#        Coordinates are 0-based with end exclusive (Python slice convention).
#
#   2) <out_prefix>.fasta
#        The same unique stretches as FASTA records, with wrapped sequence lines.
#
#   3) <out_prefix>.hits.tsv
#        Per-k-mer hits showing WHERE the query matches the reference:
#          query_id   q_start0   q_end0_excl   kmer   ref_file   ref_header   ref_pos0   strand
#        - ref_header corresponds to the FASTA header token captured by the parser
#          (currently the first whitespace-delimited token after '>').
#        - ref_pos0 is the 0-based start coordinate of the k-mer within that reference record.
#        - strand indicates whether the k-mer matched the reference forward ('+')
#          or via the reverse-complement key ('-').
#        - max_hits_per_kmer limits how many reference hits are printed per query k-mer to avoid
#          extremely large output for repetitive k-mers.
#
# Example:
#   python3 unique_stretches.py \
#     --query "TCCTGGCTCAGGACGAACGCT..." \
#     --ref_dir /path/to/reference_fastas \
#     -k 35 \
#     --min_len 50 \
#     --out_prefix my_probe_candidates
#
# Notes / caveats:
#   - The reference index is more memory-intensive than a simple set because it stores
#     provenance (file, header, position, strand) for each k-mer occurrence.
#   - If the reference database is very large, this k-mer->hit-list index may require
#     substantial RAM; consider a disk-backed approach or a Bloom filter if needed.
#   - FASTA header capture currently uses only the first token after '>' (split on whitespace).
#     If you want the entire header line preserved, change:
#         header = line[1:].split()[0]
#     to:
#         header = line[1:].strip()


import argparse
import gzip
import sys
from pathlib import Path

DNA_COMP = str.maketrans("ACGTN", "TGCAN")

def rc(seq: str) -> str:
    # reverse-complement DNA (T, not U)
    return seq.translate(DNA_COMP)[::-1]

def normalize_query(seq: str) -> str:
    # Accept RNA or DNA input; convert to DNA alphabet
    seq = seq.strip().upper().replace("U", "T")
    # keep only letters; if you want to allow gaps/spaces, strip them here
    return "".join([c for c in seq if c.isalpha()])

def open_text_maybe_gz(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")

def iter_fasta_sequences(path: Path):
    """
    Yields (header, sequence) from a FASTA file (optionally gzipped).
    Streaming parser: does not load the whole file into memory.
    """
    with open_text_maybe_gz(path) as fh:
        header = None
        chunks = []
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks).upper()
                header = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if header is not None:
            yield header, "".join(chunks).upper()

def build_reference_kmer_index(ref_dir: Path, k: int, exts=(".fa", ".fasta", ".fna", ".ffn", ".fas")) -> dict:
    """
    Walk ref_dir recursively, read FASTA files (and .gz), build an index:
        kmer -> list of hits
    Each hit is a tuple: (ref_file, header, pos0, strand)
    We store both kmer (strand '+') and rc(kmer) (strand '-') so queries match either.
    NOTE: This is more memory heavy than a plain set.
    """
    index = {}  # kmer -> list[(ref_file, header, pos0, strand)]
    files = []

    for p in ref_dir.rglob("*"):
        if not p.is_file():
            continue
        name = p.name.lower()
        base = name[:-3] if name.endswith(".gz") else name
        if any(base.endswith(e) for e in exts):
            files.append(p)

    if not files:
        raise SystemExit(f"No FASTA-like files found under: {ref_dir}")

    for idx, f in enumerate(sorted(files), 1):
        print(f"[ref {idx}/{len(files)}] indexing {f}", file=sys.stderr)
        for header, seq in iter_fasta_sequences(f):
            seq = seq.upper().replace("U", "T")
            n = len(seq)
            if n < k:
                continue
            for pos in range(0, n - k + 1):
                kmer = seq[pos:pos+k]
                if "N" in kmer or any(c not in "ACGT" for c in kmer):
                    continue

                # forward orientation hit
                index.setdefault(kmer, []).append((str(f), header, pos, "+"))

                # reverse-complement key so queries match either strand
                rk = rc(kmer)
                index.setdefault(rk, []).append((str(f), header, pos, "-"))

    print(f"[ref] total unique {k}-mers stored: {len(index):,}", file=sys.stderr)
    return index

def find_unique_stretches(query_seq: str, ref_kmers: set, k: int, min_len: int):
    """
    Compute boolean array absent[i] for each k-mer starting at i in the query:
      absent[i] = (query[i:i+k] not in ref_kmers)
    Then convert runs of True into maximal base-interval stretches:
      run i..j corresponds to bases [i, j+k) (end exclusive)
    Only keep stretches with length >= min_len.
    Returns list of (start0, end_excl0, subseq).
    """
    q = normalize_query(query_seq)
    n = len(q)
    if n < k:
        return []

    absent = [False] * (n - k + 1)
    for i in range(0, n - k + 1):
        kmer = q[i:i+k]
        if "N" in kmer or any(c not in "ACGT" for c in kmer):
            # treat ambiguous query windows as "not eligible"
            absent[i] = False
            continue
        absent[i] = (kmer not in ref_kmers)

    stretches = []
    i = 0
    while i < len(absent):
        if not absent[i]:
            i += 1
            continue
        j = i
        while j + 1 < len(absent) and absent[j + 1]:
            j += 1
        start = i
        end = j + k  # end exclusive
        subseq = q[start:end]
        if len(subseq) >= min_len:
            stretches.append((start, end, subseq))
        i = j + 1

    return stretches

def write_hits_tsv(query_seq: str, query_name: str, ref_index: dict, k: int, out_path: str, max_hits_per_kmer: int = 5):
    """
    For each query k-mer window, if present in ref_index, write where it matches:
      query_start0, query_end0_excl, kmer, ref_file, ref_header, ref_pos0, strand
    max_hits_per_kmer limits spam when a k-mer occurs many times.
    """
    q = normalize_query(query_seq)
    n = len(q)
    with open(out_path, "w") as out:
        out.write("query_id\tq_start0\tq_end0_excl\tkmer\tref_file\tref_header\tref_pos0\tstrand\n")
        if n < k:
            return
        for i in range(0, n - k + 1):
            kmer = q[i:i+k]
            if "N" in kmer or any(c not in "ACGT" for c in kmer):
                continue
            hits = ref_index.get(kmer)
            if not hits:
                continue
            # write up to max_hits_per_kmer hits
            for (rf, hdr, pos0, strand) in hits[:max_hits_per_kmer]:
                out.write(f"{query_name}\t{i}\t{i+k}\t{kmer}\t{rf}\t{hdr}\t{pos0}\t{strand}\n")


def main():
    ap = argparse.ArgumentParser(
        description="Find maximal query stretches (>=min_len) whose every k-mer window is absent from reference FASTA database."
    )
    gq = ap.add_mutually_exclusive_group(required=True)
    gq.add_argument("--query", help="Query sequence (RNA or DNA).")
    gq.add_argument("--query_fasta", help="FASTA file containing query sequence(s). Only the first record is used.")

    ap.add_argument("--ref_dir", required=True, help="Directory containing reference FASTA/FA/FNA files (can be .gz). Searched recursively.")
    ap.add_argument("-k", type=int, default=35, help="k-mer length (default: 35).")
    ap.add_argument("--min_len", type=int, default=35, help="Minimum length of output stretches (default: 35).")
    ap.add_argument("--out_prefix", default="unique_stretches", help="Output prefix (default: unique_stretches).")
    args = ap.parse_args()

    ref_dir = Path(args.ref_dir)
    if not ref_dir.exists():
        raise SystemExit(f"ref_dir not found: {ref_dir}")

    # Load query
    if args.query is not None:
        query_seq = args.query
        query_name = "query"
    else:
        qpath = Path(args.query_fasta)
        if not qpath.exists():
            raise SystemExit(f"query_fasta not found: {qpath}")
        recs = list(iter_fasta_sequences(qpath))
        if not recs:
            raise SystemExit(f"No FASTA records found in: {qpath}")
        query_name, query_seq = recs[0]

    k = args.k
    min_len = args.min_len
    if min_len < k:
        print(f"[warn] min_len < k; setting min_len = k ({k})", file=sys.stderr)
        min_len = k

    # Build reference k-mer index (kmer -> list of reference hits)
    ref_index = build_reference_kmer_index(ref_dir, k)
    ref_kmers = set(ref_index.keys())  # membership-only view for unique stretch detection


    # Find unique stretches
    stretches = find_unique_stretches(query_seq, ref_kmers, k, min_len)
    print(f"[query] length: {len(normalize_query(query_seq))} bp", file=sys.stderr)
    print(f"[query] unique stretches found: {len(stretches)}", file=sys.stderr)

    # Write outputs


    ##tsv_path = f"{args.out_prefix}.tsv"
    fa_path  = f"{args.out_prefix}.fasta"
    hits_path = f"{args.out_prefix}.hits.tsv"

    ##Path(tsv_path).parent.mkdir(parents=True, exist_ok=True)
    Path(fa_path).parent.mkdir(parents=True, exist_ok=True)
    Path(hits_path).parent.mkdir(parents=True, exist_ok=True)

    ##with open(tsv_path, "w") as tsv:
    ##    tsv.write("query_id\tstart0\tend0_excl\tlength\tsequence\n")
    ##    for idx, (s, e, seq) in enumerate(stretches, 1):
    ##        tsv.write(f"{query_name}\t{s}\t{e}\t{e-s}\t{seq}\n")

    with open(fa_path, "w") as fa:
        for idx, (s, e, seq) in enumerate(stretches, 1):
            fa.write(f">{query_name}|stretch{idx}|{s}-{e}|len={e-s}\n")
            # wrap 60 chars
            for i in range(0, len(seq), 60):
                fa.write(seq[i:i+60] + "\n")

    print(f"[out] wrote {tsv_path}", file=sys.stderr)

    # Also write a per-kmer hits table (where the query IS found in the reference)
    write_hits_tsv(query_seq, query_name, ref_index, k, hits_path, max_hits_per_kmer=5)
    print(f"[out] wrote {hits_path}", file=sys.stderr)


if __name__ == "__main__":
    main()

