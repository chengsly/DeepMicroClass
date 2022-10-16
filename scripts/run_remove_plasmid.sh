#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate tf2

python scripts/remove_plasmid.py -i data/post20_archaea -o data/clean_post20_archaea
python scripts/remove_plasmid.py -i data/post20_bacteria -o data/clean_post20_bacteria