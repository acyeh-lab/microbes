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
```install -c conda-forge ncbi-datasets-cli``` for python interface.  

Otherwise, to install binaries:
```curl -o datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/mac/datasets'
curl -o dataformat 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/mac/dataformat'
chmod +x datasets dataformat
mv datasets /usr/local/bin/
mv dataformat /usr/local/bin/
datasets version
dataformat version```
## /usr/local/bin/ is part of your system’s default $PATH environment variable, which is a colon-separated list of directories that your shell searches to find executables.

