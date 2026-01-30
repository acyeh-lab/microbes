## To search for genome hits:
**This function counts how many microbial genomes in this folder contain this DNA sequence (on either strand)**
- Must take into account files being *.fna.gz
- Also, must search for REVERSE COMPLEMENT as well! (forward genome is arbitrary)



```
pattern="atgaggattgatatattaattggacatactagtttttttcatcaaaccagtagagataacttccttcactatctcaatgaggaagaaataaaacgctatgatcagtttcattttgtgagtgataaagaactctatattttaagccgtatcctgctcaaaacagcactaaaaagatatcaacctgatgtctcattacaatcatggcaatttagtacgtgcaaatatggcaaaccatttatagtttttcctcagttggcaaaaaagattttttttaacctttcccatactatagatacagtagccgttgctattagttctcactgcgagcttggtgtcgatattgaacaaataagagatttagacaactcttatctgaatatcagtcagcatttttttactccacaggaagctactaacatagtttcacttcctcgttatgaaggtcaattacttttttggaaaatgtggacgctcaaagaagcttacatcaaatatcgaggtaaaggcctatctttaggactggattgtattgaatttcatttaacaaataaaaaactaacttcaaaatatagaggttcacctgtttatttctctcaatggaaaatatgtaactcatttctcgcattagcctctccactcatcacccctaaaataactattgagctatttcctatgcagtcccaactttatcaccacgactatcagctaattcattcgtcaaatgggcagaattga"

> hits.txt  # truncate / create output file

for f in *.fna.gz; do
  if zcat "$f" | awk -v pat="$pattern" '
    BEGIN {
      RS=">"; FS="\n"

      # Uppercase the pattern
      up = ""
      for (i = 1; i <= length(pat); i++) {
        up = up toupper(substr(pat, i, 1))
      }
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
    }

    NR > 1 {
      seq = ""
      for (i = 2; i <= NF; i++) {
        seq = seq toupper($i)
      }
      if (index(seq, pat) || index(seq, pat_rc)) {
        exit 0
      }
    }
    END { exit 1 }
  '; then
    echo "$f" >> hits.txt
  fi
done
```
