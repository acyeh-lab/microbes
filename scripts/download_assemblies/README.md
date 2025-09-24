The script "predownload_assemblies.py" script downloads all the GCF/GCA assemblies in a given file (we use a search from bv-brc for Akkermansia as as example).

The script "download_ncbi_datasets.py" script downloads all complete, annotated genomes, proteoms, and gff3 files for a given organism. The gff3 download enables 16S rRNA extraction.


## Example run of Akkermansia genomes:

### If using a .csv file of entries
```
sbatch predownload_assemblies.sh \
  /fh/fast/hill_g/Albert/Collaboration-Microbiome/blast_results/Akkermansia_muciniphila.csv\
  /fh/fast/hill_g/Albert/Collaboration-Microbiome/blast_results/Akkermansia_muciniphila\
  12 \
  /fh/fast/hill_g/Albert/Collaboration-Microbiome/blast_results/Akkermansia_muciniphila/manifest.csv
```
### If downloading directly from NCBI using command line input for species:
```
sbatch download_ncbi_dataset.sh \
  "akkermansia muciniphila" \
  /fh/fast/hill_g/Albert/Collaboration-Microbiome/NCBI/Akkermansia_muciniphila
```

### For manual download without using scripts
To download off of command line, activate environment of interest, and install ncbi-datasets-cli:
```install -c conda-forge ncbi-datasets-cli``` for python interface.  

Otherwise, to install binaries:
```
curl -o datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/mac/datasets'
curl -o dataformat 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/mac/dataformat'
chmod +x datasets dataformat
mv datasets /usr/local/bin/
mv dataformat /usr/local/bin/
datasets version
dataformat version
```
**Note that /usr/local/bin/ is part of your system’s default $PATH environment variable, which is a colon-separated list of directories that your shell searches to find executables.

Then to install genomes from command line (https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/genomes/download-genome/).  Examples below with parameters:
```datasets download genome taxon human --assembly-level complete```
```datasets download genome taxon human --annotated```
```datasets download genome taxon human --reference```
```datasets download genome taxon human --reference --include genome,rna,cds,protein,gtf```
For example, to download all e.coli genomes, proteoms, and gtf files that are annotated and with complete assembly:

```datasets download genome taxon "akkermansia muciniphila" --annotated --assembly-level complete --include genome,protein,gtf,gbff```
Note that gbff may contain more reliable than gtf for 16S rRNA extraction.


## Extracting Genomes, Proteomes, and 16S
Run script:
```sbatch unpack_single_organism_16S.sh /fh/fast/hill_g/Albert/Collaboration-Microbiome/NCBI/Akkermansia_muciniphila
```
to also parse 16S or
```sbatch unpack_single_organism.sh /fh/fast/hill_g/Albert/Collaboration-Microbiome/NCBI/Akkermansia_muciniphila
```


