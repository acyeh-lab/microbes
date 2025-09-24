#!/bin/bash
#SBATCH --job-name=download_organism
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --partition campus-new
# Usage:
# sbatch download_single_organism.sh "Akkermansia muciniphila" /path/to/output_dir

# Get input arguments
ORGANISM="$1"
OUTDIR="$2"

# Activate your micromamba environment
ENV_PREFIX="/home/ayeh/micromamba/envs/microbiome_genomics"
export PATH="$ENV_PREFIX/bin:$PATH"

# For debugging
echo "Using environment: $ENV_PREFIX"
which python
which datasets

# Create output directory if needed
mkdir -p "$OUTDIR"
mkdir -p logs

echo "🔍 Downloading all genomes for organism: $ORGANISM"
echo "📁 Output directory: $OUTDIR"

# Run the script
/home/ayeh/micromamba/envs/probe_design/bin/python download_ncbi_dataset.py --organism "$ORGANISM" --outdir "$OUTDIR"

echo "Done downloading $ORGANISM genomes!"

