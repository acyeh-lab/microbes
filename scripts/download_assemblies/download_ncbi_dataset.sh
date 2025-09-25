#!/bin/bash
#SBATCH --job-name=download_organism
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --partition campus-new

# Usage:
# sbatch download_ncbi_dataset.sh "Akkermansia muciniphila" /path/to/output_dir TRUE
# (third argument "TRUE" means include --reference flag)

# Get input arguments
ORGANISM="$1"
OUTDIR="$2"
REFERENCE="$3"

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

echo "Downloading all genomes for organism: $ORGANISM"
echo "Output directory: $OUTDIR"
echo "Reference flag: $REFERENCE"

# Build the Python command
PYTHON_CMD="$ENV_PREFIX/bin/python /fh/fast/hill_g/Albert/Collaboration-Microbiome/Scripts/download_ncbi_dataset.py --organism \"$ORGANISM\" --outdir \"$OUTDIR\""

# Add --reference flag if specified
if [ "$REFERENCE" == "TRUE" ]; then
    PYTHON_CMD="$PYTHON_CMD --reference"
fi

# Execute
eval $PYTHON_CMD

echo "Done downloading $ORGANISM genomes!"

