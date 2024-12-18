
To link sequence accession numbers (like GenBank accessions) to their taxonomy IDs, download NCBI's accession-to-taxid mapping files.
```
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
```

To annotate genomes:
```
module load prokka/1.14.5-gompi-2020b
```
```
prokka GCF_XXXXXXX.1.fna --outdir annotation_output --prefix annotated_genome
```
