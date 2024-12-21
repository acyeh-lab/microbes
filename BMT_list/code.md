
To link sequence accession numbers (like GenBank accessions) to their taxonomy IDs, download NCBI's accession-to-taxid mapping files.
```
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
```
Of note, this is a huge file - it contains all bacterial taxa and accession numbers associated with that taxa (all genes).


For bacterial genomes, we used NCBI Datasets (https://github.com/acyeh-lab/NCBIDatasets) to download all relevant taxa / accession numbers into the following directory:
```
/fh/fast/hill_g/Albert/Genomes_Proteomes/NCBI_Datasets/
```

We then run Prokka to obtain annotated datasets For example:
```
/fh/fast/hill_g/Albert/Genomes_Proteomes/NCBI_Datasets/Akkermansia_muciniphilia/ncbi_dataset/data/GCF_009731575.1/annotation_output
```

A description of Prokka output can be found here: "https://github.com/acyeh-lab/Prokka"

Using this pipeline, we have stored annotated datasets for all of the relevant taxa, for example, using the base directory "/fh/fast/hill_g/Albert/Genomes_Proteomes/NCBI_Datasets/":
```
Bacteroides_thetaiotaomicron/ncbi_dataset/data/GCF_014131755.1/annotation_output
Akkermansia_muciniphilia/ncbi_dataset/data/GCF_009731575.1/annotation_output
Enterococcus_faecalis/ncbi_dataset/data/GCF_000393015.1/annotation_output
```

This can assist in the search of genes, such as the ```buk``` gene.

To delete mid-level "ncbi_dataset.zip" files after unzipping:
```
find . -mindepth 2 -maxdepth 2 -type f -name "*.zip" -exec rm -v {} \;
```
