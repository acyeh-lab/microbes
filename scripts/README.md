## To search for genome hits:
- Must take into account files being *.fna.gz
- Also, must search for REVERSE COMPLEMENT as well! (forward genome is arbitrary)



```
pattern="ATGGCTGGAACCTGTAGTACGGACTTACGAGA"
for f in *.fna.gz; do
  if zcat "$f" | awk -v pat="$pattern" '
    BEGIN {
      RS=">"; FS="\n"

      # 1) Uppercase the pattern once
      up = ""
      for (i = 1; i <= length(pat); i++) {
        c  = substr(pat, i, 1)
        up = up toupper(c)
      }
      pat = up

      # 2) Build reverse complement of the pattern (also uppercase)
      rev = ""
      for (i = length(pat); i >= 1; i--) {
        c = substr(pat, i, 1)
        if      (c == "A") rc = "T"
        else if (c == "T") rc = "A"
        else if (c == "C") rc = "G"
        else if (c == "G") rc = "C"
        else               rc = c      # keep N or other symbols as-is
        rev = rev rc
      }
      pat_rc = rev
    }

    # Skip the empty first record caused by RS=">"
    NR > 1 {
      seq = ""
      # concatenate all sequence lines, uppercased
      for (i = 2; i <= NF; i++) {
        seq = seq toupper($i)
      }

      # Look for either forward or reverse complement
      if (index(seq, pat) > 0 || index(seq, pat_rc) > 0) {
        found = 1
        exit 0    # found something in this file → stop early
      }
    }

    END {
      if (found) exit 0
      else       exit 1
    }
  '; then
    echo "$f"
  fi
done | wc -l
```
