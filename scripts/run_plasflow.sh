#!/bin/sh

eval "$(conda shell.bash hook)"
conda activate plasflow

# for f in data/sampled_sequence_subsampled/*.fa
for f in data/sampled_sequence_subsampled_2000/*.fa
do
    filename=$(basename -- "$f")
    extension="${filename##*.}"
    filename="${filename%.*}"
    echo "Predicting $f"
    PlasFlow.py --input $f --output result_other_2000/result_plasflow/$filename
done