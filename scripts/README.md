# This folder contains a collection of random scripts I used at some point or another


## find_shared_regions.py
```
# Summary:
#   Identify "shared" regions in a chosen reference sequence that are conserved
#   across a directory of FASTA files (each file treated as one strain).
#
#   The script uses exact k-mer matching (default k=35). It:
#     1) Loads all *.fasta files in the given directory (each file = one strain).
#     2) For each strain, builds a UNION set of all k-mers observed anywhere in
#        that file's records (contigs/chromosomes).
#     3) Selects a reference record (sequence) from which all candidate k-mers are taken.
#     4) Counts, for each reference k-mer, in how many strains (files) that k-mer appears.
#     5) Keeps only k-mers that meet a "strain support" threshold:
#          - If --min-fraction is set: require ceil(min_fraction * #strains)
#          - Else if --require-at-least is set: require that many strains
#          - Else: require presence in ALL strains
#     6) Marks base positions on the reference that lie within any supported k-mer,
#        then extracts continuous regions on the reference where:
#          - Region length >= k
#          - Every k-mer window inside the region is supported ("shared")
#     7) Merges overlapping/contiguous regions and writes them to a TSV.
#
# What "shared" means here:
#   A k-mer is "shared" if it is present (exact string match) in at least the
#   required number of STRAINS (files). Presence is counted by strain, not by
#   number of contigs/records.
#
# Inputs (CLI arguments):
#   1) fasta_dir (positional): Directory containing *.fasta files (each file = one strain)
#   2) -k / --kmer: k-mer length, and minimum region length (default 35)
#   3) --min-fraction: fraction of strains required to contain the k-mer (0..1]
#   4) --require-at-least: integer # of strains required (ignored if --min-fraction is set)
#   5) --ref-file: choose a specific reference strain FASTA file (uses its first record)
#   6) --ref-name: substring match to choose a specific reference record name (if --ref-file not set)
#   7) -o / --out: output TSV path (default shared_regions.tsv)
#
# Reference selection (priority order):
#   1) --ref-file: first FASTA record in that file
#   2) --ref-name: first record whose constructed name contains the substring
#   3) default: first record in the first *.fasta file (alphabetical order)
#
# Output:
#   TSV file with a small header block plus tab-separated columns:
#     sequence   start   end
#   Where:
#     - 'sequence' is the reference substring for each shared region
#     - 'start' is 1-based on the reference
#     - 'end' is written as the reference slice end index (see code); interpret as
#       1-based inclusive end if using start=s+1 and end=e for the slice [s:e)
#
# Example:
#   python shared_regions.py genomes_dir/ -k 35 --min-fraction 0.9 --ref-file GCF_000001.fasta -o shared.tsv
#
# Notes / caveats:
#   - Only files ending with *.fasta are read (not .fa/.fna/.gz unless modified).
#   - Reverse complements are NOT considered; matches are strand-sensitive.
#   - "Strain presence" is based on k-mer existence anywhere in that strain file.
#   - Runtime can be high for large reference sequences because support counting
#     scales roughly with (#reference_kmers × #strains).
#   - Memory can be high because per-strain k-mer UNION sets are stored in RAM.
```

## geneome_search.sh
```
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
```



## Scratchwork - To search for genome hits:
**This function counts how many microbial genomes in this folder contain this DNA sequence (on either strand)**
- Must take into account files being *.fna.gz
- Also, must search for REVERSE COMPLEMENT as well! (forward genome is arbitrary)



```
pattern="G"

> hits.txt

for f in *.fna.gz; do
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
    echo "$f" >> hits.txt
  fi
done

```
