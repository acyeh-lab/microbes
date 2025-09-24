The script "predownload_assemblies.py" script downloads all the GCF/GCA assemblies in a given file (we use a search from bv-brc for Akkermansia as as example).

## Example run of Akkermansia genomes:

```
sbatch predownload_assemblies.sh \
  /fh/fast/hill_g/Albert/Collaboration-Microbiome/blast_results/Akkermansia_muciniphila.csv\
  /fh/fast/hill_g/Albert/Collaboration-Microbiome/blast_results/Akkermansia_muciniphila\
  12 \
  /fh/fast/hill_g/Albert/Collaboration-Microbiome/blast_results/Akkermansia_muciniphila/manifest.csv
```

To download off of command line, activate environment of interest, and install ncbi-datasets-cli:
```install -c conda-forge ncbi-datasets-cli```
