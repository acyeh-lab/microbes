#!/usr/bin/env bash
#SBATCH --job-name=genome_search
#SBATCH --output=logs/unique_stretches_%j.out
#SBATCH --error=logs/unique_stretches_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --partition campus-new

set -euo pipefail

# Summary:
#   SLURM wrapper to run unique_stretches.py for identifying maximal "unique"
#   query stretches (>= min_len) whose every k-mer window (default k=35) is absent
#   from a reference FASTA database. Also writes a per-kmer hits table indicating
#   where query k-mers DO match the reference (file, FASTA header, position, strand).
#
# Inputs:
#   Mode A (literal query string):
#     sbatch run_unique_stretches.sbatch --query "<SEQ>" REF_DIR [OUT_PREFIX] [K] [MIN_LEN]
#
#   Mode B (query FASTA; first record used):
#     sbatch run_unique_stretches.sbatch --query_fasta QUERY.fa REF_DIR [OUT_PREFIX] [K] [MIN_LEN]
#
# Required:
#   --query OR --query_fasta
#   REF_DIR      Directory of reference FASTA files (searched recursively; .gz allowed)
#
# Optional:
#   OUT_PREFIX   Output prefix (default: unique_stretches)
#   K            k-mer length (default: 35)
#   MIN_LEN      Minimum unique-stretch length (default: 35)
#
# Output:
#   <OUT_PREFIX>.tsv        Unique stretches (query-centric)
#   <OUT_PREFIX>.fasta     FASTA of unique stretches
#   <OUT_PREFIX>.hits.tsv  Per-kmer reference hits with FASTA headers and coordinates
#
# Example:
#   ./run_unique_stretches.sh --query "AATGGAAACAGGTGCTAATACCGCATAACAGTTTA" \
#     /fh/fast/hill_g/Albert/Collaboration-Microbiome/NCBI_data/Ref_gut_human/db/Enterococcus_faecalis/16S_rRNA \
#     my_probe_candidates 35 35
#
# Notes:
#   - Uses python3 explicitly (important on clusters where `python` is Python 2)
#   - Logs written to logs/genome_search_<jobid>.out/.err
#   - Ensure the logs/ directory exists before submission

usage() {
  cat <<'EOF'
Usage:
  sbatch run_unique_stretches.sbatch --query "<SEQ>" REF_DIR [OUT_PREFIX] [K] [MIN_LEN]
  sbatch run_unique_stretches.sbatch --query_fasta QUERY_FASTA REF_DIR [OUT_PREFIX] [K] [MIN_LEN]

Examples:
  sbatch run_unique_stretches.sbatch \
    --query "AATGGAAACAGGTGCTAATACCGCATAACAGTTTA" \
    /path/to/ref_dir \
    my_probe_candidates 35 35

  sbatch run_unique_stretches.sbatch \
    --query_fasta query.fa \
    /path/to/ref_dir \
    my_probe_candidates 40 50
EOF
}

# --- basic argument checking ---
if [[ $# -lt 3 ]]; then
  usage
  exit 2
fi

MODE="$1"
QVAL="$2"
REF_DIR="$3"
OUT_PREFIX="${4:-unique_stretches}"
K="${5:-35}"
MIN_LEN="${6:-35}"

# Ensure logs directory exists
mkdir -p logs

# Resolve script path
SCRIPT_DIR="/fh/fast/hill_g/Albert/Collaboration-Microbiome/Scripts"
PY_SCRIPT="${SCRIPT_DIR}/unique_stretches.py"

if [[ ! -f "$PY_SCRIPT" ]]; then
  echo "[ERROR] unique_stretches.py not found at: $PY_SCRIPT" >&2
  exit 2
fi

if [[ ! -d "$REF_DIR" ]]; then
  echo "[ERROR] REF_DIR is not a directory: $REF_DIR" >&2
  exit 2
fi

# Build base command
CMD=(python3 "$PY_SCRIPT"
     --ref_dir "$REF_DIR"
     -k "$K"
     --min_len "$MIN_LEN"
     --out_prefix "$OUT_PREFIX")

case "$MODE" in
  --query)
    CMD+=(--query "$QVAL")
    ;;
  --query_fasta)
    if [[ ! -f "$QVAL" ]]; then
      echo "[ERROR] query_fasta not found: $QVAL" >&2
      exit 2
    fi
    CMD+=(--query_fasta "$QVAL")
    ;;
  *)
    echo "[ERROR] First argument must be --query or --query_fasta" >&2
    usage
    exit 2
    ;;
esac

echo "[INFO] Job started on $(hostname) at $(date)" >&2
echo "[INFO] Running command:" >&2
printf ' %q' "${CMD[@]}" >&2
echo >&2

"${CMD[@]}"

echo "[INFO] Job finished at $(date)" >&2
echo "[INFO] Outputs:" >&2
echo "  ${OUT_PREFIX}.tsv" >&2
echo "  ${OUT_PREFIX}.fasta" >&2
echo "  ${OUT_PREFIX}.hits.tsv" >&2

