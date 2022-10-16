#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate rfplasmid

for f in $HOME/project/DeepMicrobeFinder/data/sampled_sequence_subsampled/*.fa
do
    filename=$(basename -- "$f")
    extension="${filename##*.}"
    filename="${filename%.*}"
    echo "Predicting $f"
    rfplasmid --input data/rfplasmid/${filename} --threads 16 --out $HOME/project/DeepMicrobeFinder/data/result_rfplasmid/${filename} --species Generic
    # mkdir data/rfplasmid/${filename}
    # cp $f data/rfplasmid/${filename}
done