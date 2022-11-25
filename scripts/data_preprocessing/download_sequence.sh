#!/bin/bash

working_dir=/home/tianqi/project/DeepMicrobeFinder/data

cd $working_dir

if [ ! -f "genome_updater.sh" ]; then
    wget --quiet --show-progress https://raw.githubusercontent.com/pirovc/genome_updater/master/genome_updater.sh
    chmod +x genome_updater.sh
fi

############################################################################################################
# Download bacteria, archaea, and eukaryota genomes from NCBI RefSeq and GenBank
############################################################################################################

./genome_updater.sh -d "refseq,genbank" -l "complete genome" -f "genomic.fna.gz" -t 12 -o "sequence/pre20_archaea" -T '2157' -E 20191230 -R 10 
./genome_updater.sh -d "refseq,genbank" -l "complete genome" -f "genomic.fna.gz" -t 12 -o "sequence/post20_archaea" -T '2157' -D 20200101 -R 10 
./genome_updater.sh -d "refseq" -l "complete genome" -f "genomic.fna.gz" -t 12 -o "sequence/pre20_bacteria" -T '2' -E 20191230 -R 10 
./genome_updater.sh -d "refseq" -l "complete genome" -f "genomic.fna.gz" -t 12 -o "sequence/post20_bacteria" -T '2' -D 20200101 -R 10 

./genome_updater.sh -d "refseq,genbank" -f "genomic.fna.gz" -t 12 -o "sequence/pre20_eukaryote" -E 20191230 -T '554915,554296,3027,33682,33634,33630,543769,260810,5752,2611341,2763,3041,38254' -R 10 -l "complete genome"

./genome_updater.sh -d "refseq,genbank" -f "genomic.fna.gz" -t 12 -o "sequence/post20_eukaryote" -D 20200101 -T '554915,554296,3027,33682,33634,33630,543769,260810,5752,2611341,2763,3041,38254' -R 10 -l "complete genome"

############################################################################################################
# Download virus-host db
############################################################################################################
vhdb_dir=/home/tianqi/dataset