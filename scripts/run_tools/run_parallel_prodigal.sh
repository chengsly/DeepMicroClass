#!/bin/sh

eval "$(conda shell.bash hook)"
conda activate tiara

for f in data/sampled_sequence_subsampled/*.fa
do
    filename=$(basename -- "$f")
    extension="${filename##*.}"
    filename="${filename%.*}"
    echo "Predicting $f"
    scripts/parallel_prodigal.sh $f "data/result_whokaryote/${filename}" 16 data/tmp
done