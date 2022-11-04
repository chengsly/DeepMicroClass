#!/bin/sh

eval "$(conda shell.bash hook)"
conda activate tiara

# for f in data/sampled_sequence_subsampled/*.fa
for f in data/sampled_sequence_subsampled_2000/*.fa
do
    filename=$(basename -- "$f")
    extension="${filename##*.}"
    filename="${filename%.*}"
    echo "Predicting $f"
    tiara -i $f -o "data/result_tiara_2000/${filename}.txt" -t 16 -m 0 # -p 0.2 0.3
done