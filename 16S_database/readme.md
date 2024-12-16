# Constructing a 16S rRNA database
This branch describes how to build a microbial 16S database for blast searching and MSA

# Downloading NCBI's Microbial 16S rRNA DB
NCBI has a project called "Microbial 16S RNA sequences from the RefSeq Targeted Loci project": "https://www.ncbi.nlm.nih.gov/books/NBK62345/table/blast_ftp_site.T.preformatted_refseq_n",  
In command line / shell, run the following:
```
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/16S_ribosomal_RNA.tar.gz
```
In the same directory, then run:
```
tar -xzvf 16S_ribosomal_RNA.tar.gz
```



