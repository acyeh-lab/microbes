#!/usr/bin/env python3
import argparse, csv, os, re, sys, shutil, time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import subprocess

ACC_RE = re.compile(r"\bGC[AF]_\d+\.\d+\b")

def try_run(cmd, cwd=None):
    p = subprocess.run(
        cmd, cwd=str(cwd) if cwd else None,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    return p.returncode == 0, p.stdout, p.stderr, p.returncode

# --- NEW: robust extractor that ignores headers and just finds GCA/GCF anywhere ---
def extract_accessions(table: Path, acc_col: str | None) -> list[str]:
    """
    Read the file as raw text (CSV/TSV/plain), strip BOM, and extract all
    GCA_/GCF_ tokens via regex. No header required.
    """
    text = table.read_text(encoding="utf-8-sig", errors="ignore")
    # handle non-breaking/zero-width chars if present
    text = text.replace("\u00A0", " ").replace("\u200b", "")
    accs = sorted(set(re.findall(r"\bGC[AF]_\d+\.\d+\b", text)))
    return accs

def concat_fnas_to(fna_files: list[Path], out_fna: Path) -> None:
    out_fna.parent.mkdir(parents=True, exist_ok=True)
    with out_fna.open("w") as w:
        for p in sorted(fna_files):
            with p.open("r") as r:
                shutil.copyfileobj(r, w)

# --- NEW: make the datasets call tolerant of CLI differences ---
def datasets_fetch(acc_dir: Path, acc: str) -> tuple[bool, str]:
    """
    Try to download using 'datasets download genome accession' with a couple of
    --include variants. Returns (ok, message).
    """
    # prefer genome+gbff when supported; fall back to genome only
    for include in ("genome,gbff", "genome"):
        ok, out, err, _ = try_run([
            "datasets", "download", "genome", "accession", acc,
            "--no-progressbar", "--include", include
        ], cwd=acc_dir)
        if ok:
            return True, ""
        # if the error complains about invalid argument for --include, try next option
        combined = (err or "") + (out or "")
        if "invalid argument" in combined and "include" in combined:
            continue
        # any other error: return immediately
        return False, (err.strip() or out.strip() or "datasets failed")
    # exhausted options
    return False, "datasets: could not find a supported --include option"

def download_one(acc: str, cache_dir: Path, retries: int = 2) -> dict:
    """
    Returns dict: {accession, status, fasta, size_bytes, message}
    """
    acc_dir = cache_dir / acc
    out_fna = acc_dir / f"{acc}_genomic.fna"
    if out_fna.exists() and out_fna.stat().st_size > 0:
        return {
            "accession": acc, "status": "cached", "fasta": str(out_fna),
            "size_bytes": out_fna.stat().st_size, "message": ""
        }

    acc_dir.mkdir(parents=True, exist_ok=True)
    attempt = 0
    last_msg = ""
    while attempt <= retries:
        attempt += 1

        # 1) download (zip name is auto by datasets)  --- CHANGED to use datasets_fetch ---
        ok, msg = datasets_fetch(acc_dir, acc)
        if not ok:
            last_msg = f"datasets failed: {msg}"
            time.sleep(1.0 * attempt)
            continue

        # 2) find newest zip
        zips = sorted(acc_dir.glob("*.zip"), key=lambda p: p.stat().st_mtime, reverse=True)
        if not zips:
            last_msg = "no zip found after datasets download"
            time.sleep(1.0 * attempt)
            continue

        # 3) unzip (overwrite if re-run)
        ok, out, err, _ = try_run(["unzip", "-o", str(zips[0]), "-d", str(acc_dir)])
        if not ok:
            last_msg = f"unzip failed: {err.strip() or out.strip()}"
            time.sleep(1.0 * attempt)
            continue

        # 4) collect genomic .fna(s)
        data_root = acc_dir / "ncbi_dataset" / "data"
        fna_files = list(data_root.rglob("*_genomic.fna"))
        if not fna_files:
            # fallback to any .fna
            fna_files = list(data_root.rglob("*.fna"))
        if not fna_files:
            last_msg = "no *_genomic.fna (or .fna) files found in datasets output"
            time.sleep(1.0 * attempt)
            continue

        try:
            concat_fnas_to(fna_files, out_fna)
        except Exception as e:
            last_msg = f"concat failed: {e}"
            time.sleep(1.0 * attempt)
            continue

        if out_fna.exists() and out_fna.stat().st_size > 0:
            return {
                "accession": acc, "status": "ok", "fasta": str(out_fna),
                "size_bytes": out_fna.stat().st_size, "message": ""
            }
        else:
            last_msg = "empty output fasta after concat"

    return {
        "accession": acc, "status": "error", "fasta": str(out_fna),
        "size_bytes": out_fna.stat().st_size if out_fna.exists() else 0,
        "message": last_msg
    }

def main():
    ap = argparse.ArgumentParser(description="Pre-download all assemblies (GCA_/GCF_) into a local cache for later BLAST use.")
    ap.add_argument("--table", required=True, help="Input CSV/TSV with Assembly Accession column or entries")
    ap.add_argument("--acc-col", default=None, help="Column name to read (default: autodetect).")
    ap.add_argument("--cache-dir", default="genome_cache", help="Where to store downloads and FASTAs.")
    ap.add_argument("--workers", type=int, default=6, help="Parallel downloads.")
    ap.add_argument("--manifest", default="manifest.csv", help="Write a summary manifest here.")
    ap.add_argument("--debug", action="store_true", help="Print debug info about parsed accessions.")
    args = ap.parse_args()

    table = Path(args.table)
    cache_dir = Path(args.cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)

    accs = extract_accessions(table, args.acc_col)
    if not accs:
        print("No GCA_/GCF_ accessions found.", file=sys.stderr)
        sys.exit(2)

    if args.debug:
        print(f"DEBUG: Found {len(accs)} accessions. First 5: {accs[:5]}", flush=True)

    print(f"Found {len(accs)} assemblies. Downloading with {args.workers} workers…")

    results = []
    with ThreadPoolExecutor(max_workers=max(1, args.workers)) as ex:
        fut2acc = {ex.submit(download_one, acc, cache_dir): acc for acc in accs}
        done = 0
        for fut in as_completed(fut2acc):
            acc = fut2acc[fut]
            try:
                res = fut.result()  # {accession, status, fasta, size_bytes, message}
            except Exception as e:
                res = {"accession": acc, "status": "error", "fasta": "", "size_bytes": 0, "message": str(e)}
            results.append(res)
            done += 1
            if (done % 5 == 0) or (done == len(accs)):
                print(f"Progress: {done}/{len(accs)}", flush=True)
            # Print a per-accession status line for Slurm logs
            if res["status"] in ("ok", "cached"):
                print(f"[{done}/{len(accs)}] {res['status'].upper()}  {res['accession']}  ->  {res['fasta']}  ({res['size_bytes']} bytes)", flush=True)
            else:
                print(f"[{done}/{len(accs)}] ERROR  {res['accession']}  :  {res['message']}", flush=True)

    # write manifest
    man_path = Path(args.manifest)
    with man_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["accession","status","fasta","size_bytes","message"])
        w.writeheader()
        w.writerows(results)

    # summary
    ok = sum(1 for r in results if r["status"] in ("ok","cached"))
    err = sum(1 for r in results if r["status"] == "error")
    print(f"Done. OK/Cached: {ok}   Errors: {err}   Manifest: {man_path}")

if __name__ == "__main__":
    main()
