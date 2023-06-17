#!/bin/bash

# This script is used to run the whole process of data preprocessing.

# Activate conda environment pytorch113
eval "$(conda shell.bash hook)"
conda activate pytorch113

# The working directory is the directory where the data is stored.
working_dir=/home/tianqi/project/DeepMicrobeFinder/data
cd $working_dir
./download_sequence.sh

python remove_plasmid.py -i pre20_archaea/pre20_archaea.fa -o clean_pre20/pre20_archaea.fa
python remove_plasmid.py -i pre20_bacteria/pre20_bacteria.fa -o clean_pre20/pre20_bacteria.fa
python remove_plasmid.py -i pre20_eukaryote/pre20_euk.fa -o clean_pre20/pre20_euk.fa

python remove_plasmid.py -i post20_archaea/post20_archaea.fa -o clean_post20/post20_archaea.fa
python remove_plasmid.py -i post20_bacteria/post20_bacteria.fa -o clean_post20/post20_bacteria.fa
python remove_plasmid.py -i post20_eukaryote/post20_euk.fa -o clean_post20/post20_euk.fa

./calculate_mash.sh

mkdir single_fasta
python split_fasta.py -i clean_pre20/pre20_archaea.fa -o single_fasta
python split_fasta.py -i clean_pre20/pre20_bacteria.fa -o single_fasta
python split_fasta.py -i clean_pre20/pre20_euk.fa -o single_fasta
python split_fasta.py -i clean_post20/post20_archaea.fa -o single_fasta
python split_fasta.py -i clean_post20/post20_bacteria.fa -o single_fasta
python split_fasta.py -i clean_post20/post20_euk.fa -o single_fasta