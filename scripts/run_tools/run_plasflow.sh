#!/bin/sh

eval "$(conda shell.bash hook)"
conda activate plasflow

# for f in data/sampled_sequence_subsampled/*.fa
for f in data/sampled_sequence/10000/*.fa
do
    filename=$(basename -- "$f")
    extension="${filename##*.}"
    filename="${filename%.*}"
    echo "Predicting $f"
    PlasFlow.py --input $f --output data/result_others/result_plasflow/$filename
done