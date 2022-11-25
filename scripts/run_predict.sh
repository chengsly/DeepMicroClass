#!/bin/sh

eval "$(conda shell.bash hook)"
conda activate tf2

for f in data/sampled_sequence/*.fa
do
    echo "Predicting $f"
    python predict.py -i $f -e one-hot -d models/one-hot-models -m hybrid
done


