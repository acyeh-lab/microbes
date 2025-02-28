# Constructing a 16S rRNA database
This branch describes how to build a microbial 16S database for blast searching and MSA.
Note I used this DB below to provide to 10x to design Xenium probes

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
ml BLAST+/2.14.0-gompi-2022b 
blastn -query query.fasta -db 16S_ribosomal_RNA -out results.txt
```

## Find all species listed
Directory where dataase is stored: /fh/fast/hill_g/Albert/Bacterial_Taxonomy/16S_rRNA_DB
```
blastdbcmd -db 16S_ribosomal_RNA -entry all -outfmt "%a %T %S" > species_summary.txt
```
For this file, we shoudl get 27,159 entries.
%a is accession number (file)
%T is official taxonomy ID
%S is species name
%d is date of last update (helpful to find the most recent one if multiple hits) - if this info is not available, output will be in ASN format.

To get the number of unique taxa, can do:
```
cut -d' ' -f2 species_summary.txt | sort | uniq | wc -l
```
We get 21,116 here. Not so bad (out of 27,159 entries).


## Search for a species
```
grep "Bacteroides thetaiotaomicron" species_summary.txt
```
Sample output:
```
NR_112944.1 818 Bacteroides thetaiotaomicron
NR_112142.1 818 Bacteroides thetaiotaomicron
NR_074277.1 818 Bacteroides thetaiotaomicron
```
Note, these are 3 separate accession numbers for the same taxa (818) - difficult to tell which one is more recent as the %d info is not available for this blast db.  So we will keep all of them.

## Get the 16S sequence for a taxon:
```
blastdbcmd -db 16S_ribosomal_RNA -entry NR_112944.1 -outfmt "%f" > NR_112944.1.fasta
```
Note that "NR" prefix stands for non-coding RNA sequence used in refSeq annotation.
```
grep "Bifidobacterium" species_summary.txt | grep -o "NR_[^ ]*" | paste -sd "," - | sed 's/,/, /g'
```
This returns the list of all "NR_" accession numbers for the genus Bifidobacterium

## Comprehensive Human Gut Reference
https://gmrepo.humangut.info/home

Deposited here: /fh/fast/hill_g/Albert/Bacterial_Taxonomy/Human_Gut_Reference

