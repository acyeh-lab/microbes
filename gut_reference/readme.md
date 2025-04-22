## Human gut reference database
- Compiled from: https://gmrepo.humangut.info/home
- Deposited here (generated 2/3/25): ```/fh/fast/hill_g/Albert/Bacterial_Taxonomy/Human_Gut_Reference/GMREPO_relative_abundance_of_all_species_genus_in_all_phenotypes_summary.tsv```
- The file "human_gut_reference_database.csv" is a curated dataset from the above reference.  It includes all bacteria that occurred in at least 2 samples of a given disease group and only taxa that occurred in > 0.01%.
- The file "taxon_ids.txt" contains just the taxon ids, and is used to download the corresponding genomes below.

## Download references using NCBI taxon ID
- Use NCBI Datasets command line tools (CLI): https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/
Download datasets, download dataformat, and make them executable:
```
curl -o datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets'
curl -o dataformat 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/dataformat'
chmod +x datasets dataformat
```

- Depoisted here (generated 2/19/25):
- ```/fh/fast/hill_g/Albert/Bacterial_Taxonomy/Human_Gut_Reference```
- Run command (after getting command line tools to work - need to create local ~/bin folder as need to move the datasets and dataformat binaries to a directory that is in your PATH, so that you can run them from anywhere without specifying the full path).
```
mkdir -p ~/bin
mv datasets dataformat ~/bin/
chmod +x ~/bin/datasets ~/bin/dataformat
echo $PATH
echo 'export PATH=$HOME/bin:$PATH' >> ~/.bashrc
source ~/.bashrc
```
- Then run the following script to download all of the taxon IDs to the local directory
```
while read taxid; do
    if [ -d "genomes/$taxid" ]; then
        echo "Skipping Taxon ID: $taxid (already downloaded)"
        continue
    fi
    
    echo "Downloading genome for Taxon ID: $taxid"
    datasets download genome taxon "$taxid" --dehydrated --include genome --include cds --include seq-report --filename "$taxid.zip"
    
    if [ -f "$taxid.zip" ]; then
        unzip -o "$taxid.zip" -d "genomes/$taxid"
    else
        echo "Error: Failed to download genome for Taxon ID: $taxid"
    fi
done < taxon_ids.txt
```
- Move all IDs into the folder "genome". (Make the genome folder first)
- Then rehydrate all datasets:
```
datasets rehydrate --directory genomes/
```
