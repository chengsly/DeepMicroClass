#!/bin/bash

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate genomad

# Run genomad
# for f in data/sampled_sequence_subsampled_2000/*.fa
for f in data/sampled_sequence_subsampled_10000/*.fa
do
    filename=$(basename -- "$f")
    extension="${filename##*.}"
    filename="${filename%.*}"
    echo "Predicting $f"
    genomad end-to-end -t 32 --enable-score-calibration --cleanup $f $HOME/project/DeepMicrobeFinder/data/result_others_10000/result_genomad/ /home/tianqi/software/genomad/genomad_db
done