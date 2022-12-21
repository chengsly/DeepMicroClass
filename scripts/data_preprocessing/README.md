# Data preprocessing
This directory contains the scripts used to download and preprocess the data used in this project.
## Data downloading
The dataset we used to train and evaluate the project consists of five classes: 

- Prokaryote
    - Archaea (RefSeq)
    - Bacteria (RefSeq)
- Eukaryote
    - RefSeq
    - MMETSP
- Plasmid
    - PLSDB
- Virus
    - Virus-Host DB
        - Prokaryotic virus
        - Eukaryotic virus
    - IMG/VR

The Archaea, Bacteria were downloaded with [genome_updater](https://github.com/pirovc/genome_updater).
Archaea were downloaded from RefSeq and Genbank, and bacteria were downloaded only from RefSeq or the dataset would be too large.
Microbial eukaryotes were downloaded from RefSeq, and additional marine eukaryotic host re-assemblies were downloaded from the [MMETSP](https://zenodo.org/record/1212585) project.
Plasmid sequences and metadata were downloaded from [PLSDB](https://ccb-microbe.cs.uni-saarland.de/plsdb).
Virus sequences and metadata were downloaded from [Virus-Host DB](https://www.genome.jp/ftp/db/virushostdb/).
In addition to the experimentally verified virus genomes, we also retrived the IMG/VR virome dataset.

The whole downloading scripts can be found in `download_sequence.sh`.

## Plasmid cleaning
The bacteria and archaea genomes can have plasmid sequences, to remove the plasmid from the genomes, run `python remove_plasmid.py`.
The script first checks the sequence description, and removes any sequence with "plasmid" in it, then the sequence ids of bacteria and archaea are cross compared with the sequence ids of plasmids in PLSDB, and any sequences with overlapping ids were removed.

## Mash distance calculation
To avoid sequences too similar to each other in training and testing set, we use [Mash](https://github.com/marbl/Mash/releases) to calculate the Mash distances between each pairs of sequences within a particular class.

## Splitting training, validation and test set

For sequences from NCBI, `download_sequence.sh` downloads the archaea, bacteria and eukaryote sequences separated by submission date 01/01/2020, and stored into different directories.
For PLSDB, the plasmid sequences are separated by submission date 01/01/2020 according to the metadata, and output to two separate fasta file, the corresponding script is in `process_plsdb.py`.

## Converting sequences into onehot
