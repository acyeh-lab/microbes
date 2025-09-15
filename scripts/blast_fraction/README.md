This function returns what fraction of strains an input sequence is found in from "query.fa'. This can be used to find how good a probe coverage is. The output is named the same as the .fa file, which should contain 1 probe.  To query mutltiple probes, have to submit job with different ".fa' file.

```
sbatch blast_fraction_by_assembly.sh \
  --table /fh/fast/hill_g/Albert/Collaboration-Microbiome/blast_results/Akkermansia_muciniphila.csv \
  --sequence /fh/fast/hill_g/Albert/Collaboration-Microbiome/blast_results/akk_muc1.fa \
  --cache-dir /fh/fast/hill_g/Albert/Collaboration-Microbiome/blast_results/Akkermansia_muciniphila \
  --out /fh/fast/hill_g/Albert/Collaboration-Microbiome/blast_results/Akkermansia_muciniphila/results.csv \
  --min-pident 100 --min-qcov 100 --max-evalue 1e-5 \
  --env-prefix /home/ayeh/micromamba/envs/microbiome_genomics
```

## BLAST Filtering Thresholds

These flags control when an alignment **counts as a hit** in the results.  
A hit is recorded only if it passes **all three** thresholds.

### `--min-pident <percent>`
Require **≥ N% sequence identity** across the aligned region.

- **Definition:** `pident = (number_of_identical_positions / alignment_length) × 100`
- **Example:** If 900 of 1000 aligned bases match, `pident = 90%`.
- **Typical ranges:**  
  - Same species/strain (DNA): `90–98`  
  - Cross-strain proteins (tblastn): `35–60`  
  - Short amplicons (e.g., 16S regions): `97–100`

---

### `--min-qcov <percent>`
Require **≥ N% query coverage**, i.e., aligned length as a percentage of the full query length.

- **Definition:** `qcov = (alignment_length / query_length) × 100`
- **Example:** Query length = 1000 bp, aligned length = 820 bp → `qcov = 82%` ✅.
- **Why it matters:** Prevents tiny partial matches from counting.

---

### `--max-evalue <value>`
Require **E-value ≤ X** (expected number of chance alignments this good).

- **Lower = more stringent**. E-value is sensitive to query length and database size.
- **Typical values:** `1e-3` (lenient) → `1e-10` (stringent).  
  Defaults here: `1e-5`.

---

## Recommended Settings

Pick a preset based on your use case:

| Use case                             | Suggested flags                                         |
|-------------------------------------|---------------------------------------------------------|
| Same/near strain (DNA)              | `--min-pident 95 --min-qcov 90 --max-evalue 1e-10`     |
| Cross-strain protein (tblastn)      | `--min-pident 40 --min-qcov 70 --max-evalue 1e-5`      |
| Short amplicon (e.g., 200–400 bp)   | `--min-pident 97 --min-qcov 95 --max-evalue 1e-3`      |
| Sensitive, allow distant matches    | `--min-pident 80 --min-qcov 70 --max-evalue 1e-3`      |
| **Exact full-length match only**    | `--min-pident 100 --min-qcov 100 --max-evalue 1e-20`   |

> **Note:** For protein queries, this pipeline uses `tblastn` against nucleotide genomes by default. Identity thresholds for proteins can be lower than for DNA.

---

## Examples

**DNA query (default stringent settings):**
```bash
python blast_fraction_by_assembly.py \
  --table Akkermansia_muciniphila.csv \
  --sequence query.fa \
  --cache-dir genome_cache \
  --skip-download \
  --out results.csv \
  --min-pident 90 --min-qcov 80 --max-evalue 1e-5
