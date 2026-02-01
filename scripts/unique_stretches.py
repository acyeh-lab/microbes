#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Comments / Implementation Notes
#
# 1) Memory safety (large reference collections)
#    This script intentionally DOES NOT build a full reference k-mer index for the
#    entire database. Instead, it:
#      - extracts k-mers only from the QUERY (<= len(query)-k+1 kmers)
#      - streams reference FASTA files once
#      - records hits only for those query k-mers (bounded by --max_hits_per_kmer)
#    This avoids MemoryError when scanning large genome collections
#    (e.g. >10,000 bacterial assemblies).
#
# 2) Strand handling
#    Reference scanning checks both:
#      - forward k-mer in the reference ("+")
#      - reverse-complement matches via a precomputed mapping ("-")
#    Query sequences are NOT reverse-complemented.
#    Strand in the hits output reflects how the match was observed in the reference.
#
# 3) Definition of "unique"
#    A stretch of the query is considered UNIQUE if every k-mer window fully inside
#    that stretch is ABSENT from the reference database.
#    Practically:
#      - we mark each query k-mer window as present/absent
#      - we merge consecutive absent windows into maximal base-interval stretches
#    Unique stretches are reported only if their length >= --min_len (and --min_len is
#    forced to be >= k to avoid degenerate cases).
#
# 4) Coordinate conventions
#    - All coordinates are 0-based
#    - End positions are end-exclusive
#    - Query hit windows are [q_start0, q_end0_excl) where q_end0_excl = q_start0 + k
#    This matches Python slicing semantics and avoids off-by-one ambiguity.
#
# 5) FASTA headers
#    The full FASTA header line (everything after '>') is preserved and reported
#    in the hits output. This is important for:
#      - multi-copy genes (e.g. multiple 16S rRNA loci)
#      - distinguishing contigs / replicons in multi-FASTA reference files
#
# 6) Ambiguous bases
#    Any k-mer containing:
#      - 'N'
#      - non-ACGT IUPAC ambiguity codes (R,Y,S,W,K,M,…)
#    is skipped in both query and reference to avoid false matches.
#
# 7) Output files
#    This script produces ONLY:
#      - <out_prefix>.fasta      (unique stretches)
#      - <out_prefix>.hits.tsv   (where query k-mers match the reference)
#    No TSV summary of unique stretches is written by design.
#
# 8) Output size controls
#    Repetitive k-mers can match many times across genomes/contigs. To keep the
#    hits output bounded, hits are capped per k-mer by --max_hits_per_kmer.
#
# 9) Parallelism + progress reporting
#    Reference scanning is parallelized across files using ProcessPoolExecutor.
#    Implementation detail:
#      - one worker task is submitted per file
#      - the parent process merges results and prints a single progress bar line as
#        each file completes
#    This provides true per-file progress without log spam from worker processes.
#
# 10) Reference file discovery (non-recursive by default)
#    list_reference_files() currently scans ONLY files directly under --ref_dir
#    (uses Path.iterdir), not subdirectories. If you want recursive discovery,
#    switch back to Path.rglob("*") and filter by extension.
#
# 11) Intended use cases
#    - Probe / primer design
#    - Identifying species- or strain-specific regions
#    - Screening candidate sequences against large microbial databases
#
# -----------------------------------------------------------------------------

import argparse
import gzip
import sys
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

DNA_COMP = str.maketrans("ACGTN", "TGCAN")

def rc(seq: str) -> str:
    return seq.translate(DNA_COMP)[::-1]

def normalize_query(seq: str) -> str:
    seq = seq.strip().upper().replace("U", "T")
    return "".join([c for c in seq if c.isalpha()])

def render_progress(done: int, total: int, width: int = 30) -> str:
    if total <= 0:
        return "[?]"
    frac = done / total
    filled = int(round(width * frac))
    bar = "#" * filled + "-" * (width - filled)
    return f"[{bar}] {done}/{total} ({frac*100:5.1f}%)"

def open_text_maybe_gz(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")

def iter_fasta_sequences(path: Path):
    """
    Yields (header, sequence) from a FASTA file (optionally gzipped).
    Header is the full line after '>' (no truncation).
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
                header = line[1:].strip()  # keep full header line
                chunks = []
            else:
                chunks.append(line)
        if header is not None:
            yield header, "".join(chunks).upper()

def query_kmers_set(query_seq: str, k: int) -> set:
    q = normalize_query(query_seq)
    out = set()
    if len(q) < k:
        return out
    for i in range(0, len(q) - k + 1):
        km = q[i:i+k]
        if "N" in km or any(c not in "ACGT" for c in km):
            continue
        out.add(km)
    return out

def build_rc_to_query(query_kmers: set) -> dict:
    # key: ref_kmer (which equals rc(query_kmer)) -> value: query_kmer
    return {rc(qk): qk for qk in query_kmers}


def list_reference_files(ref_dir: Path, exts=(".fa", ".fasta", ".fna", ".ffn", ".fas")) -> list[Path]:
    files = []
    for p in ref_dir.iterdir():   # <-- NOT recursive
        if not p.is_file():
            continue
        name = p.name.lower()
        base = name[:-3] if name.endswith(".gz") else name
        if any(base.endswith(e) for e in exts):
            files.append(p)

    if not files:
        raise SystemExit(f"No FASTA-like files found directly under: {ref_dir}")

    return sorted(files)

def scan_one_file_for_query_kmers(
    fstr: str,
    query_kmers: set,
    rc_to_query: dict,
    k: int,
    max_hits_per_kmer: int,
):
    """
    Worker: scan a single FASTA(/.gz) file, returning:
      present: set of query kmers found in this file
      hits: dict query_kmer -> list[(ref_file, header, pos0, strand)] (capped)
    """
    present = set()
    hits = {}  # only store kmers that actually hit in this file

    f = Path(fstr)
    for header, seq in iter_fasta_sequences(f):
        seq = seq.upper().replace("U", "T")
        n = len(seq)
        if n < k:
            continue

        for pos in range(0, n - k + 1):
            kmer = seq[pos:pos+k]
            if "N" in kmer or any(c not in "ACGT" for c in kmer):
                continue

            # forward match: ref kmer == query kmer
            if kmer in query_kmers:
                present.add(kmer)
                lst = hits.setdefault(kmer, [])
                if len(lst) < max_hits_per_kmer:
                    lst.append((fstr, header, pos, "+"))

            # reverse-strand match:
            # if reference kmer equals rc(query_kmer), report as strand '-'
            qk = rc_to_query.get(kmer)
            if qk is not None:
                present.add(qk)
                lst = hits.setdefault(qk, [])
                if len(lst) < max_hits_per_kmer:
                    lst.append((fstr, header, pos, "-"))

    return present, hits

def find_unique_stretches(query_seq: str, present_query_kmers: set, k: int, min_len: int):
    """
    Mark each query k-mer window as absent if it is NOT in present_query_kmers.
    Convert runs of absent windows into maximal base-interval stretches.
    Returns list of (start0, end0_excl, subseq).
    """
    q = normalize_query(query_seq)
    n = len(q)
    if n < k:
        return []

    absent = [False] * (n - k + 1)
    for i in range(0, n - k + 1):
        kmer = q[i:i+k]
        if "N" in kmer or any(c not in "ACGT" for c in kmer):
            absent[i] = False
            continue
        absent[i] = (kmer not in present_query_kmers)

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
        end = j + k  # end-exclusive interval in base coordinates
        subseq = q[start:end]
        if len(subseq) >= min_len:
            stretches.append((start, end, subseq))
        i = j + 1

    return stretches

def write_unique_fasta(stretches, query_name: str, out_prefix: str):
    fa_path = f"{out_prefix}.fasta"
    Path(fa_path).parent.mkdir(parents=True, exist_ok=True)

    with open(fa_path, "w") as fa:
        for idx, (s, e, seq) in enumerate(stretches, 1):
            fa.write(f">{query_name}|stretch{idx}|{s}-{e}|len={e-s}\n")
            for i in range(0, len(seq), 60):
                fa.write(seq[i:i+60] + "\n")

    return fa_path


def write_hits_from_hitsdict(query_seq: str, query_name: str, hits_by_kmer: dict, k: int, out_path: str):
    q = normalize_query(query_seq)
    n = len(q)
    out_dir = Path(out_path).parent
    out_dir.mkdir(parents=True, exist_ok=True)

    with open(out_path, "w") as out:
        out.write("query_id\tq_start0\tq_end0_excl\tkmer\tref_file\tref_header\tref_pos0\tstrand\n")
        if n < k:
            return
        for i in range(0, n - k + 1):
            kmer = q[i:i+k]
            if "N" in kmer or any(c not in "ACGT" for c in kmer):
                continue
            hits = hits_by_kmer.get(kmer, [])
            for (rf, hdr, pos0, strand) in hits:
                out.write(f"{query_name}\t{i}\t{i+k}\t{kmer}\t{rf}\t{hdr}\t{pos0}\t{strand}\n")

def parallel_scan_reference_for_query_kmers(
    ref_dir: Path,
    query_kmers: set,
    k: int,
    workers: int = 1,
    max_hits_per_kmer: int = 5,
    progress: bool = True,
):
    files = list_reference_files(ref_dir)
    rc_to_query = build_rc_to_query(query_kmers)

    total = len(files)
    print(f"[ref] files: {total}  workers: {workers}", file=sys.stderr)

    # single-process path
    if workers <= 1 or total <= 1:
        present_all = set()
        hits_all = {km: [] for km in query_kmers}
        done = 0

        for f in files:
            present, hits = scan_one_file_for_query_kmers(
                str(f), query_kmers, rc_to_query, k, max_hits_per_kmer
            )
            present_all |= present
            for km, lst in hits.items():
                space = max_hits_per_kmer - len(hits_all[km])
                if space > 0:
                    hits_all[km].extend(lst[:space])

            done += 1
            if progress:
                print(f"\r[ref] {render_progress(done, total)}  last={f.name}", end="", file=sys.stderr, flush=True)

        if progress:
            print("", file=sys.stderr)
        return present_all, hits_all

    # multi-process path: submit ONE future per file (true per-file progress)
    present_all = set()
    hits_all = {km: [] for km in query_kmers}

    done = 0
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {
            ex.submit(
                scan_one_file_for_query_kmers,
                str(f), query_kmers, rc_to_query, k, max_hits_per_kmer
            ): f
            for f in files
        }

        for fut in as_completed(futs):
            f = futs[fut]
            present, hits = fut.result()

            present_all |= present
            for km, lst in hits.items():
                space = max_hits_per_kmer - len(hits_all[km])
                if space > 0:
                    hits_all[km].extend(lst[:space])

            done += 1
            if progress:
                print(f"\r[ref] {render_progress(done, total)}  last={f.name}", end="", file=sys.stderr, flush=True)

    if progress:
        print("", file=sys.stderr)

    return present_all, hits_all

def main():
    ap = argparse.ArgumentParser(
        description="Find maximal query stretches (>=min_len) whose every k-mer window is absent from reference FASTA database. Also report where query k-mers hit in the reference."
    )
    gq = ap.add_mutually_exclusive_group(required=True)
    gq.add_argument("--query", help="Query sequence (RNA or DNA).")
    gq.add_argument("--query_fasta", help="FASTA file containing query sequence(s). Only the first record is used.")

    ap.add_argument("--ref_dir", required=True, help="Directory containing reference FASTA/FA/FNA files (can be .gz).")
    ap.add_argument("-k", type=int, default=35, help="k-mer length (default: 35).")
    ap.add_argument("--min_len", type=int, default=35, help="Minimum length of output stretches (default: 35).")
    ap.add_argument("--out_prefix", default="unique_stretches", help="Output prefix (default: unique_stretches).")
    ap.add_argument("--max_hits_per_kmer", type=int, default=5, help="Cap hits written per k-mer (default: 5).")
    ap.add_argument("--workers", type=int, default=1, help="Number of processes for reference scanning (default: 1).")
    ap.add_argument("--no_progress", action="store_true", help="Disable progress bar output.")


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

    q_kmers = query_kmers_set(query_seq, k)
    if not q_kmers:
        raise SystemExit(f"Query is shorter than k={k} or contains no valid A/C/G/T {k}-mers.")

    present_kmers, hits_by_kmer = parallel_scan_reference_for_query_kmers(
        ref_dir=ref_dir,
        query_kmers=q_kmers,
        k=k,
        workers=args.workers,
        max_hits_per_kmer=args.max_hits_per_kmer,
        progress=(not args.no_progress),
    )

    stretches = find_unique_stretches(query_seq, present_kmers, k, min_len)

    qlen = len(normalize_query(query_seq))
    print(f"[query] length: {qlen} bp", file=sys.stderr)
    print(f"[query] valid {k}-mers: {len(q_kmers)}", file=sys.stderr)
    print(f"[ref] query-kmers present in reference: {len(present_kmers)} / {len(q_kmers)}", file=sys.stderr)
    print(f"[query] unique stretches found: {len(stretches)}", file=sys.stderr)

    fa_path = write_unique_fasta(stretches, query_name, args.out_prefix)
    print(f"[out] wrote {fa_path}", file=sys.stderr)

    hits_path = f"{args.out_prefix}.hits.tsv"
    write_hits_from_hitsdict(query_seq, query_name, hits_by_kmer, k, hits_path)
    print(f"[out] wrote {hits_path}", file=sys.stderr)

if __name__ == "__main__":
    main()
