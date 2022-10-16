#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate vs2

for f in data/sampled_sequence_subsampled_2000/*.fa
do
    filename=$(basename -- "$f")
    extension="${filename##*.}"
    filename="${filename%.*}"
    echo "Predicting $f"
    virsorter run -w result_other_2000/result_virsorter2/${filename} -i $f -j 16 --provirus-off --max-orf-per-seq 10 --include-groups dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae all
done