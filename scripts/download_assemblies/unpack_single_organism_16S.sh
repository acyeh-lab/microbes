#!/bin/bash
#SBATCH --job-name=unpack_ncbi
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --partition campus-new
# Usage:
# sbatch unpack_single_organism.sh /path/to/downloaded_dataset_dir
# ALSO INCLUDES 16S rRNA

# Get input argument
DATASET_DIR="$1"

# Path to environment and script
ENV_PREFIX="/home/ayeh/micromamba/envs/microbiome_genomics"
SCRIPT_PATH="/fh/fast/hill_g/Albert/Collaboration-Microbiome/Scripts/unpack_ncbi_dataset_16S.py"

# Update PATH to use micromamba environment
export PATH="$ENV_PREFIX/bin:$PATH"

# Create logs directory if it doesn't exist
mkdir -p logs

# For debugging/logging
echo "Environment: $ENV_PREFIX"
which python

echo "Changing to dataset directory: $DATASET_DIR"
cd "$DATASET_DIR" || { echo "ERROR: Failed to change directory to $DATASET_DIR"; exit 1; }

echo "Running unpack script on: $DATASET_DIR"

# Run the script
"$ENV_PREFIX/bin/python" "$SCRIPT_PATH"

echo "Finished unpacking $DATASET_DIR"

