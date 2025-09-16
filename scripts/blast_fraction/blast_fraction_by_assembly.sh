#!/usr/bin/env bash
# Run blast_fraction_by_assembly.py against pre-downloaded genomes.
#
#SBATCH --job-name=blast_fraction
#SBATCH --partition=campus-new
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail

# -------- defaults (override with flags below) --------
TABLE=""
SEQUENCE=""
CACHE_DIR="genome_cache"
OUT="results.csv"   # if left as default/blank, we auto-name from --sequence + thresholds
MIN_PIDENT="90"
MIN_QCOV="80"
MAX_EVALUE="1e-5"
ENV_PREFIX="/home/ayeh/micromamba/envs/microbiome_genomics"
SCRIPT="/fh/fast/hill_g/Albert/Collaboration-Microbiome/blast_fraction_by_assembly.py"

usage() {
  cat <<USAGE
Required:
  --table PATH           CSV/TSV/plain text with GCA_/GCF_
  --sequence PATH        Query FASTA (DNA -> blastn; protein -> tblastn)

Optional:
  --cache-dir DIR        Directory containing <ACC>/<ACC>_genomic.fna (default: $CACHE_DIR)
  --out PATH             Results CSV (default: <cache-dir>/<sequence-stem>_p<minp>_q<minq>.csv)
  --min-pident N         % identity threshold (default: $MIN_PIDENT)
  --min-qcov N           % query coverage threshold (default: $MIN_QCOV)
  --max-evalue X         e-value threshold (default: $MAX_EVALUE)
  --env-prefix DIR       Conda/mamba env prefix (default: $ENV_PREFIX)
  --script PATH          Path to blast_fraction_by_assembly.py (default: $SCRIPT)
  -h|--help              Show this help

Tip: You can override Slurm resources at submit time, e.g.
  sbatch -p campus-new -c 16 --mem=32G -t 12:00:00 blast_fraction_by_assembly.sbatch.sh ...
USAGE
}

# -------- parse flags --------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --table) TABLE="$2"; shift 2;;
    --sequence) SEQUENCE="$2"; shift 2;;
    --cache-dir) CACHE_DIR="$2"; shift 2;;
    --out) OUT="$2"; shift 2;;
    --min-pident) MIN_PIDENT="$2"; shift 2;;
    --min-qcov) MIN_QCOV="$2"; shift 2;;
    --max-evalue) MAX_EVALUE="$2"; shift 2;;
    --env-prefix) ENV_PREFIX="$2"; shift 2;;
    --script) SCRIPT="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    --) shift; break;;
    *) echo "Unknown arg: $1"; usage; exit 2;;
  esac
done

mkdir -p logs
mkdir -p "$CACHE_DIR"

# required args
[[ -n "$TABLE" && -n "$SEQUENCE" ]] || { echo "ERROR: --table and --sequence are required"; usage; exit 2; }

# helper to format percentages nicely (95.0 -> 95 ; 97.50 -> 97.5)
fmt_pct() {
  local x="$1"
  printf "%s" "$x" | sed -E 's/^([0-9]+)\.0+$/\1/; s/^([0-9]+\.[0-9]*[1-9])0+$/\1/; s/^([0-9]+)\.$/\1/'
}

# -------- auto-name OUT from --sequence if user didn't override --------
if [[ "$OUT" == "results.csv" || -z "$OUT" ]]; then
  seq_base="$(basename "$SEQUENCE")"
  seq_stem="${seq_base%.*}"  # drop extension (.fa/.fasta/.faa)
  p_sfx="$(fmt_pct "$MIN_PIDENT")"
  q_sfx="$(fmt_pct "$MIN_QCOV")"
  OUT="${CACHE_DIR%/}/${seq_stem}_p${p_sfx}_q${q_sfx}.csv"
fi
mkdir -p "$(dirname "$OUT")"

echo "TABLE      : $TABLE"
echo "SEQUENCE   : $SEQUENCE"
echo "CACHE_DIR  : $CACHE_DIR"
echo "OUT        : $OUT"
echo "MIN_PIDENT : $MIN_PIDENT"
echo "MIN_QCOV   : $MIN_QCOV"
echo "MAX_EVALUE : $MAX_EVALUE"
echo "ENV_PREFIX : $ENV_PREFIX"
echo "SCRIPT     : $SCRIPT"

# -------- use env directly (no activation needed) --------
export PATH="$ENV_PREFIX/bin:$PATH"
export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# sanity checks
[[ -x "$ENV_PREFIX/bin/python" ]] || { echo "Missing $ENV_PREFIX/bin/python"; exit 127; }
command -v makeblastdb >/dev/null 2>&1 || { echo "BLAST+ (makeblastdb) not found in $ENV_PREFIX/bin"; exit 127; }
[[ -f "$TABLE" ]] || { echo "Missing --table: $TABLE"; exit 2; }
[[ -f "$SEQUENCE" ]] || { echo "Missing --sequence: $SEQUENCE"; exit 2; }
[[ -d "$CACHE_DIR" ]] || { echo "Missing --cache-dir: $CACHE_DIR"; exit 2; }
[[ -f "$SCRIPT" ]] || { echo "Missing --script: $SCRIPT"; exit 2; }

# run
"$ENV_PREFIX/bin/python" -u "$SCRIPT" \
  --table "$TABLE" \
  --sequence "$SEQUENCE" \
  --cache-dir "$CACHE_DIR" \
  --skip-download \
  --out "$OUT" \
  --min-pident "$MIN_PIDENT" \
  --min-qcov "$MIN_QCOV" \
  --max-evalue "$MAX_EVALUE"

echo "Done. Results: $OUT"

