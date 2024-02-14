#!/bin/sh

eval "$(conda shell.bash hook)"
conda activate pytorch

for f in data/sampled_sequence/10000/*.fa
# for f in data/sampled_sequence_subsampled_2000/*.fa
# for f in /home/shengwei/Viromes/04_read_mapping_per_sample/assemblies_per_sample/*.fasta
# for f in /home/tianqi/dataset/IMGVR_v2/sources/*.fasta
do
    echo "Predicting $f"
    # python predict.py -i $f -e one-hot -d models/one-hot-models -m hybrid
    # python predict.py -i $f -e one-hot -d models/one-hot-models -m single -l 2000
    # python predict_pytorch.py -i $f -e one-hot -d data/new_model -m hybrid -o data/virome/result_dmf
    # python predict_pytorch.py -i $f -e one-hot -d data/new_model -m hybrid -o data/virome/IMGVR
    # python predict_pytorch.py -i $f -e one-hot -d /home/tianqi/project/DeepMicrobeFinder/model.ckpt -m hybrid
    python predict_pytorch.py -i $f -e one-hot -d /home/tianqi/project/DeepMicrobeFinder/data/pt_logs/log_2023-10-28-16_05/checkpoint/epoch=299-step=153600-val_f1=0.923-val_acc=0.923.ckpt -m hybrid
done

python scripts/parse_result/parse_dmc.py