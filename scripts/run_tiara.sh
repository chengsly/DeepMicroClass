#!/bin/sh

eval "$(conda shell.bash hook)"
conda activate tiara

for f in data/sampled_sequence/*.fa
do
    filename=$(basename -- "$f")
    extension="${filename##*.}"
    filename="${filename%.*}"
    echo "Predicting $f"
    tiara -i $f -o "data/result_tiara/${filename}.txt" -t 16
done