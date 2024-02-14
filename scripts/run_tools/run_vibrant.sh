#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate vibrant

# for f in $HOME/project/DeepMicrobeFinder/data/sampled_sequence_subsampled/*.fa
for f in $HOME/project/DeepMicrobeFinder/data/sampled_sequence/10000/*.fa
# for f in /home/tianqi/dataset/IMGVR_v2/sources/*.fasta
do
    filename=$(basename -- "$f")
    extension="${filename##*.}"
    filename="${filename%.*}"
    echo "Predicting $f"
    python $HOME/software/VIBRANT/VIBRANT_run.py -i $f -t 16 -folder $HOME/project/DeepMicrobeFinder/data/result_others/result_vibrant/
    # python $HOME/software/VIBRANT/VIBRANT_run.py -i $f -t 16 -folder $HOME/project/DeepMicrobeFinder/data/virome/vibrant
done