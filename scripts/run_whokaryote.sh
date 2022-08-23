#!/bin/sh

eval "$(conda shell.bash hook)"
conda activate tiara

for f in data/sampled_sequence/*.fa
do
    filename=$(basename -- "$f")
    extension="${filename##*.}"
    filename="${filename%.*}"
    echo "Predicting $f"
    whokaryote.py --contigs $f --outdir "data/result_whokaryote/${filename}" --threads 16
done