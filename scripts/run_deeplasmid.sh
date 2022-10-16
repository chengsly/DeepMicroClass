#!/bin/bash

for f in data/sampled_sequence_subsampled/*.fa
do
    filename=$(basename -- "$f")
    extension="${filename##*.}"
    filename="${filename%.*}"
    echo "Predicting $f"
    mkdir -p data/result_deeplasmid/$filename
    docker run -it -v $f:/srv/jgi-ml/classifier/dl/in.fasta -v $HOME/project/DeepMicrobeFinder/data/result_deeplasmid/$filename:/srv/jgi-ml/classifier/dl/outdir billandreo/deeplasmid feature_DL_plasmid_predict.sh in.fasta outdir
done