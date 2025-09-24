## Sam Minot has develoeped pipeline to identify conserved sequences using whole-genome data

1) Download from NCBI genome browser (https://www.ncbi.nlm.nih.gov/datasets/genome/)
2) Filter out to only include genomes that are annotated and have complete assembly.
3) Run ```python3 unpack_ncbi_dataset.py``` which unpacks the genomes and proteomes in the directory that contains "ncbi_dataset" download.  Or can run shell script "unpack_ncbi_dataset.sh" (I modified this to run from shell):
  ```sbatch unpack_single_organism.sh /fh/fast/hill_g/Albert/Collaboration-Microbiome/NCBI/Collinsella_aerofaciens``` 
5) Then upload data into Cirro after going into the genomes and proteomes folder: ```CIRRO_BASE_URL=fredhutch.cirro.bio cirro upload -i```
6) After done uploaded, click on the genome dataset and run "Build Pangenome".
7) Now run "find-conserved-pangenome markers" after building the pangenome:
8) ```git clone https://github.com/FredHutch/find-conserved-pangenome-markers```
9) Load the uv module then type ```bash edit.sh```.  Note that I ran this on local enviroment (my laptop).  The compute itself uses AWS nodes.
10) Go through the prompts, log-in to the hutch (don't need to be VPN), select your datasets.
11) Variables include "minimal proportion of genomes" - set to 1 (this means that the gene has to occur in at least 100% of genoems). Minimum major allele frequency - set to 0.9 (means that the sequence has to be at least 90% match).  
