#!/usr/bin/env python3
import argparse, csv, os, re, sys, subprocess, shutil
from pathlib import Path
from typing import List, Optional, Tuple

def is_dna(seq: str) -> bool:
    import re
    seq = re.sub(r"\s+", "", seq).upper()
    dna_set = set("ACGTNRYKMSWBDHVU")
    bad = sum(1 for c in seq if c not in dna_set)
    return (bad / max(1, len(seq))) <= 0.05

def read_sequence(seq_path: Path) -> str:
    if not seq_path.exists():
        raise FileNotFoundError(f"Sequence file not found: {seq_path}")
    txt = seq_path.read_text().splitlines()
    if any(line.startswith(">") for line in txt):
        return "".join(l.strip() for l in txt if not l.startswith(">"))
    return "".join(l.strip() for l in txt)

def run(cmd: List[str], cwd: Optional[Path] = None) -> None:
    p = subprocess.run(cmd, cwd=str(cwd) if cwd else None,
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"Command failed: {' '.join(cmd)}\nSTDOUT:\n{p.stdout}\nSTDERR:\n{p.stderr}")

def try_run(cmd: List[str], cwd: Optional[Path] = None) -> Tuple[bool,str,str,int]:
    p = subprocess.run(cmd, cwd=str(cwd) if cwd else None,
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    return (p.returncode == 0, p.stdout, p.stderr, p.returncode)

# ---- Extraction tuned for your sheet ----
def extract_accessions(table_path: Path, acc_col_arg: Optional[str]) -> List[str]:
    # auto-detect delimiter by extension; default to csv
    delim = "," if table_path.suffix.lower() != ".tsv" else "\t"
    with table_path.open(newline="") as f:
        reader = csv.DictReader(f, delimiter=delim)
        header = reader.fieldnames or []
        if not header:
            raise ValueError("No header row detected in the input table.")

        # Prefer explicit arg; else prefer 'Assembly Accession'
        preferred = [
            "Assembly Accession", "assembly accession",
            "assembly_accession", "assembly", "accession"
        ]
        if acc_col_arg:
            acc_col = acc_col_arg
        else:
            norm = {h.strip().lower(): h for h in header}
            acc_col = None
            for p in preferred:
                if p.lower() in norm:
                    acc_col = norm[p.lower()]
                    break

        regex = re.compile(r"\bGC[AF]_\d+\.\d+\b")
        accessions = set()

        if acc_col:
            for row in reader:
                m = regex.findall((row.get(acc_col, "") or ""))
                accessions.update(m)
        else:
            # fallback: scan all columns for GCA_/GCF_
            rows = list(reader)
            for row in rows:
                for v in row.values():
                    accessions.update(regex.findall(str(v)))

        return sorted(accessions)

def ensure_genome_fasta_for_accession(acc: str, cache_dir: Path) -> Path:
    acc_dir = cache_dir / acc
    fasta_out = acc_dir / f"{acc}_genomic.fna"
    if fasta_out.exists() and fasta_out.stat().st_size > 0:
        return fasta_out
    acc_dir.mkdir(parents=True, exist_ok=True)
    ok, out, err, _ = try_run([
        "datasets","download","genome","accession", acc,
        "--no-progressbar","--include","genome,genbank"
    ], cwd=acc_dir)
    if not ok:
        raise RuntimeError(f"datasets download failed for {acc}:\n{out}\n{err}")
    # unzip newest zip
    import glob
    zips = sorted(glob.glob(str(acc_dir / "*.zip")), key=lambda p: Path(p).stat().st_mtime, reverse=True)
    if not zips:
        raise RuntimeError(f"No ZIP produced for {acc}.")
    ok, out, err, _ = try_run(["unzip","-o", zips[0], "-d", str(acc_dir)])
    if not ok:
        raise RuntimeError(f"unzip failed for {acc}:\n{out}\n{err}")
    # concat *_genomic.fna
    from pathlib import Path as P
    fna_candidates = list((acc_dir / "ncbi_dataset" / "data").rglob("*_genomic.fna"))
    if not fna_candidates:
        fna_candidates = list((acc_dir / "ncbi_dataset" / "data").rglob("*.fna"))
    if not fna_candidates:
        raise RuntimeError(f"No genomic FASTA found for {acc}.")
    with fasta_out.open("w") as w:
        for p in sorted(fna_candidates):
            w.write(P(p).read_text())
    return fasta_out

def ensure_blast_db(fasta: Path, dbtype: str) -> Path:
    prefix = fasta.with_suffix("")
    needed = ["nhr","nin","nsq"] if dbtype=="nucl" else ["phr","pin","psq"]
    if all((prefix.with_suffix("." + e)).exists() for e in needed):
        return prefix
    run(["makeblastdb","-in",str(fasta),"-dbtype",dbtype,"-out",str(prefix)])
    return prefix

def has_hit(query_fasta: Path, db_prefix: Path, program: str,
            min_pident: float, min_qcov: float, max_evalue: float) -> bool:
    outfmt = "6 sacc pident length qlen qstart qend sstart send evalue bitscore"
    ok, out, err, _ = try_run([
        program,"-query",str(query_fasta),"-db",str(db_prefix),
        "-outfmt",outfmt,"-evalue",str(max_evalue),
        "-max_target_seqs","5","-num_threads",str(os.cpu_count() or 1)
    ])
    if not ok:
        return False
    for line in out.strip().splitlines():
        parts = line.split("\t")
        if len(parts) < 10: 
            continue
        _, pident, length, qlen, *_rest, evalue, _bits = parts
        try:
            pident = float(pident); length = float(length); qlen = float(qlen); evalue = float(evalue)
        except ValueError:
            continue
        qcov = (length/qlen)*100.0 if qlen>0 else 0.0
        if pident >= min_pident and qcov >= min_qcov and evalue <= max_evalue:
            return True
    return False

def write_temp_query_fasta(seq: str, tmp_dir: Path) -> Path:
    tmp_dir.mkdir(parents=True, exist_ok=True)
    q = tmp_dir / "query.fasta"
    with q.open("w") as f:
        f.write(">query\n")
        for i in range(0, len(seq), 80):
            f.write(seq[i:i+80] + "\n")
    return q

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", required=True, help="CSV/TSV with an 'Assembly Accession' column")
    ap.add_argument("--sequence", required=True, help="Query sequence (FASTA or raw)")
    ap.add_argument("--out", default="blast_accession_results.csv")
    ap.add_argument("--acc-col", default=None, help="Override column name if needed")
    ap.add_argument("--cache-dir", default="genome_cache")
    ap.add_argument("--min-pident", type=float, default=90.0)
    ap.add_argument("--min-qcov", type=float, default=80.0)
    ap.add_argument("--max-evalue", type=float, default=1e-5)
    ap.add_argument("--skip-download", action="store_true")
    args = ap.parse_args()

    table = Path(args.table)
    seq_path = Path(args.sequence)
    cache_dir = Path(args.cache_dir); cache_dir.mkdir(parents=True, exist_ok=True)

    seq = read_sequence(seq_path)
    dna = is_dna(seq)
    program = "blastn" if dna else "tblastn"  # use tblastn for protein vs nucleotide dbs
    dbtype = "nucl"

    query_fasta = write_temp_query_fasta(seq, cache_dir / "_tmp")

    accs = extract_accessions(table, args.acc_col)
    if not accs:
        print("No GCA_/GCF_ assembly accessions found.", file=sys.stderr)
        sys.exit(2)

    print(f"Found {len(accs)} assembly accessions.")
    results = []
    hits = 0
    for i, acc in enumerate(accs, 1):
        status = "fail"; msg = ""
        try:
            if args.skip_download:  # allow manual pre-population of FASTAs
                fasta = cache_dir / acc / f"{acc}_genomic.fna"
                if not fasta.exists():
                    raise FileNotFoundError(f"Missing cached FASTA: {fasta}")
            else:
                fasta = ensure_genome_fasta_for_accession(acc, cache_dir)
            db_prefix = ensure_blast_db(fasta, dbtype)
            hit = has_hit(query_fasta, db_prefix, program,
                          args.min_pident, args.min_qcov, args.max_evalue)
            status = "hit" if hit else "no_hit"
            hits += int(hit)
        except Exception as e:
            status = "error"; msg = str(e)[:500]
        results.append({"assembly_accession": acc, "status": status, "error": msg})
        if i % 5 == 0 or i == len(accs):
            print(f"Progress: {i}/{len(accs)} processed")

    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["assembly_accession","status","error"])
        w.writeheader(); w.writerows(results)

    frac = hits / len(accs)
    print(f"Done. Hits: {hits}/{len(accs)} (fraction = {frac:.3f})")
    print(f"Wrote: {args.out}")

if __name__ == "__main__":
    main()

