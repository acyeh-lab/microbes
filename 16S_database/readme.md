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

# Understanding NCBI BLAST Database Files
The **16S_ribosomal_RNA.tar.gz** file from NCBI contains a BLAST database for 16S ribosomal RNA sequences. After extracting the archive, you'll find various files with specific suffixes, each serving a unique purpose in the BLAST database structure.
## File Suffixes and Their Roles

| **Suffix**      | **File Name Example**        | **Description**                                                                 |
|------------------|------------------------------|---------------------------------------------------------------------------------|
| `.ndb`          | `16S_ribosomal_RNA.ndb`      | Metadata about the BLAST database, used internally by BLAST.                   |
| `.nhr`          | `16S_ribosomal_RNA.nhr`      | Stores headers (e.g., sequence IDs, descriptions) for the sequences.           |
| `.nin`          | `16S_ribosomal_RNA.nin`      | Index file for locating sequence information in the `.nsq` file.               |
| `.nnd`          | `16S_ribosomal_RNA.nnd`      | Neighbor index data for advanced BLAST searches.                               |
| `.nni`          | `16S_ribosomal_RNA.nni`      | Complementary indexing details for neighbor data, used with `.nnd`.            |
| `.nog`          | `16S_ribosomal_RNA.nog`      | Offset group file for efficient database access.                               |
| `.nos`          | `16S_ribosomal_RNA.nos`      | Offset source file for optimizing BLAST searches.                              |
| `.not`          | `16S_ribosomal_RNA.not`      | Offset target file, used alongside `.nos`.                                     |
| `.nsq`          | `16S_ribosomal_RNA.nsq`      | Contains the actual nucleotide sequence data in the database.                  |
| `.ntf`          | `16S_ribosomal_RNA.ntf`      | Stores taxonomy information for the sequences, enabling taxonomic queries.     |
| `.nto`          | `16S_ribosomal_RNA.nto`      | Taxonomy offsets, used with `.ntf` to store taxonomic data efficiently.        |
| `.tar` or `.tar.gz` | `16S_ribosomal_RNA.tar.gz` | Archive containing all the database files, downloaded from NCBI.               |

## How to Use the Database in BLAST
These files form a BLAST database that can be queried using BLAST tools. For example:

```
blastn -query query.fasta -db 16S_ribosomal_RNA -out results.txt
```

# Find all species listed
```
blastdbcmd -db 16S_ribosomal_RNA -entry all -outfmt "%a %T %S" > species_summary.txt
```
