#!/bin/sh

eval "$(conda shell.bash hook)"
conda activate pytorch

# for f in data/sampled_sequence/*.fa
# for f in data/sampled_sequence_subsampled_2000/*.fa
# for f in /home/shengwei/Viromes/04_read_mapping_per_sample/assemblies_per_sample/*.fasta
for f in /home/tianqi/dataset/IMGVR_v2/sources/*.fasta
do
    echo "Predicting $f"
    # python predict.py -i $f -e one-hot -d models/one-hot-models -m hybrid
    # python predict.py -i $f -e one-hot -d models/one-hot-models -m single -l 2000
    # python predict_pytorch.py -i $f -e one-hot -d data/new_model -m hybrid -o data/virome/result_dmf
    python predict_pytorch.py -i $f -e one-hot -d data/new_model -m hybrid -o data/virome/IMGVR
done


