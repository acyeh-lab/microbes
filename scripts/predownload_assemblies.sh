#!/usr/bin/env bash
#
# SBATCH submission script for pre-downloading assemblies with NCBI Datasets.

#SBATCH --job-name=predl_assemblies
#SBATCH --partition=campus-new
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00

# NOTE: Slurm opens these log files at submit time.
# Make sure the 'logs' directory exists *before* you run `sbatch`.
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
# If you prefer absolute paths for logs, use (uncomment & edit):
# #SBATCH --output=/fh/fast/hill_g/Albert/Bacterial_Taxonomy/logs/%x-%j.out
# #SBATCH --error=/fh/fast/hill_g/Albert/Bacterial_Taxonomy/logs/%x-%j.err

set -euo pipefail

# ------- Inputs (positional with sensible defaults) -------
TABLE=${1:-/mnt/data/Akkermansia_muciniphila.csv}
CACHE_DIR=${2:-genome_cache}
WORKERS=${3:-${SLURM_CPUS_PER_TASK:-8}}
MANIFEST=${4:-${CACHE_DIR}/manifest.csv}

# ------- Use your micromamba env directly (no activation needed) -------
ENV_PREFIX="/home/ayeh/micromamba/envs/microbiome_genomics"
export PATH="${ENV_PREFIX}/bin:${PATH}"
export PYTHONUNBUFFERED=1   # stream prints to Slurm immediately

# Create runtime dirs (this happens after job starts; logs/ must exist before sbatch)
mkdir -p "${CACHE_DIR}"

echo "TABLE      : ${TABLE}"
echo "CACHE_DIR  : ${CACHE_DIR}"
echo "WORKERS    : ${WORKERS}"
echo "MANIFEST   : ${MANIFEST}"
echo "ENV_PREFIX : ${ENV_PREFIX}"

# ------- Sanity checks -------
[ -x "${ENV_PREFIX}/bin/python" ] || { echo "Python not found at ${ENV_PREFIX}/bin/python"; exit 127; }
command -v datasets >/dev/null 2>&1 || { echo "ERROR: 'datasets' CLI not found in ${ENV_PREFIX}/bin"; exit 127; }
command -v unzip >/dev/null 2>&1 || { echo "ERROR: 'unzip' not found in ${ENV_PREFIX}/bin"; exit 127; }

# ------- Run downloader -------
"${ENV_PREFIX}/bin/python" -u /fh/fast/hill_g/Albert/Collaboration-Microbiome/predownload_assemblies.py \
  --table "${TABLE}" \
  --cache-dir "${CACHE_DIR}" \
  --workers "${WORKERS}" \
  --manifest "${MANIFEST}"

echo "Done: manifest at ${MANIFEST}"
