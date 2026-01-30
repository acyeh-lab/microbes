#!/usr/bin/env bash
#
# genome_search.sh
#
# Summary:
#   Scan a directory of gzipped FASTA genome files (*.fna.gz) for the presence
#   of a DNA sequence pattern on either the forward strand or its reverse
#   complement.
#
#   For each genome file, the script:
#     - Streams and parses FASTA records without loading entire files into memory
#     - Searches for the pattern (case-insensitive) and its reverse complement
#     - Stops scanning a genome as soon as a match is found (efficient)
#
# Inputs (via sbatch):
#   1) Directory containing *.fna.gz genome files
#   2) DNA sequence pattern to search (e.g. "AGAG")
#   3) Output filename (written inside the target directory)
#
# Output:
#   A text file where:
#     - Line 1 is the search pattern used
#     - Each subsequent line is the name of a genome file containing ≥1 match
#
# Example:
#   sbatch genome_search.sh /fh/fast/hill_g/Albert "AGAG" AGAG_hits.txt
#
# Notes:
#   - FASTA strand orientation is arbitrary; both strands are searched
#   - Script is safe for large genomes and large collections of files
#   - Designed for SLURM environments
#
#SBATCH --job-name=genome_search
#SBATCH --output=genome_search_%j.out
#SBATCH --error=genome_search_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

set -euo pipefail

# -------------------------
# Argument checking
# -------------------------
if [[ $# -ne 3 ]]; then
  echo "Usage: sbatch $0 /path/to/folder \"PATTERN\" output_filename.txt"
  echo "Example: sbatch $0 /fh/fast/hill_g/Albert \"AGAG\" AGAG_hits.txt"
  exit 2
fi

DIR="$1"
pattern="$2"
OUT="$3"

if [[ ! -d "$DIR" ]]; then
  echo "ERROR: Directory not found: $DIR"
  exit 1
fi

cd "$DIR"

# -------------------------
# Safe globbing
# -------------------------
shopt -s nullglob

# -------------------------
# Initialize output file (header = pattern)
# -------------------------
echo "$pattern" > "$OUT"

files=( *.fna.gz )
if [[ ${#files[@]} -eq 0 ]]; then
  echo "No .fna.gz files found in: $DIR"
  exit 0
fi

# -------------------------
# Main search loop
# -------------------------
for f in "${files[@]}"; do
  if zcat "$f" | awk -v pat="$pattern" '
    BEGIN {
      RS=">"; FS="\n"

      # Uppercase pattern
      up = ""
      for (i = 1; i <= length(pat); i++)
        up = up toupper(substr(pat, i, 1))
      pat = up

      # Reverse complement
      rev = ""
      for (i = length(pat); i >= 1; i--) {
        c = substr(pat, i, 1)
        if      (c == "A") rc = "T"
        else if (c == "T") rc = "A"
        else if (c == "C") rc = "G"
        else if (c == "G") rc = "C"
        else rc = c
        rev = rev rc
      }
      pat_rc = rev
      found = 0
    }

    NR > 1 {
      seq = ""
      for (i = 2; i <= NF; i++)
        seq = seq toupper($i)

      if (index(seq, pat) || index(seq, pat_rc)) {
        found = 1
        exit
      }
    }

    END { exit(found ? 0 : 1) }
  '; then
    echo "$f" >> "$OUT"
  fi
done

# -------------------------
# Summary
# -------------------------
echo "Done."
echo "Directory: $DIR"
echo "Pattern:   $pattern"
echo "Output:    $DIR/$OUT"
echo -n "Hit count (excluding header): "
tail -n +2 "$OUT" | wc -l
