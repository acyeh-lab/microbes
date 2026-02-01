#!/usr/bin/env bash
#SBATCH --job-name=genome_search
#SBATCH --output=logs/unique_stretches_%j.out
#SBATCH --error=logs/unique_stretches_%j.err
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G
#SBATCH --partition campus-new

set -euo pipefail
module purge
module load Python/3.11.5-GCCcore-13.2.0

echo "[INFO] Python: $(which python3)" >&2
python3 -V >&2

usage() {
  cat <<'EOF'
Usage:
  sbatch run_unique_stretches.sbatch --query "<SEQ>" REF_DIR [OUT_PREFIX] [K] [MIN_LEN]
  sbatch run_unique_stretches.sbatch --query_fasta QUERY_FASTA REF_DIR [OUT_PREFIX] [K] [MIN_LEN]
EOF
}

# ---- minimal arg sanity check ----
if [[ $# -lt 2 ]]; then
  usage
  exit 2
fi

MODE="$1"
QVAL="$2"
REF_DIR="${3:-}"
OUT_PREFIX="${4:-unique_stretches}"
K="${5:-35}"
MIN_LEN="${6:-35}"

# Ensure required positional REF_DIR exists
if [[ -z "${REF_DIR}" ]]; then
  echo "[ERROR] REF_DIR is required" >&2
  usage
  exit 2
fi

# Ensure logs + outputs directories exist
mkdir -p logs outputs

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

OUT_PREFIX="outputs/${OUT_PREFIX}"

CMD=(python3 "$PY_SCRIPT"
     --ref_dir "$REF_DIR"
     -k "$K"
     --min_len "$MIN_LEN"
     --out_prefix "$OUT_PREFIX"
     --workers "${SLURM_CPUS_PER_TASK:-1}")

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
echo "[INFO] Using ${SLURM_CPUS_PER_TASK:-1} worker(s)" >&2
echo "[INFO] Running command:" >&2
printf ' %q' "${CMD[@]}" >&2
echo >&2

"${CMD[@]}"

echo "[INFO] Job finished at $(date)" >&2
echo "[INFO] Outputs written:" >&2
echo "  ${OUT_PREFIX}.fasta" >&2
echo "  ${OUT_PREFIX}.hits.tsv" >&2


