#!/bin/sh

eval "$(conda shell.bash hook)"
conda activate dvf
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CONDA_PREFIX/lib/

for f in data/sampled_sequence_subsampled_2000/*.fa
# for f in /home/shengwei/Viromes/04_read_mapping_per_sample/assemblies_per_sample/*.fasta
do
    filename=$(basename -- "$f")
    extension="${filename##*.}"
    filename="${filename%.*}"
    echo "Predicting $f"
    python $HOME/software/DeepVirFinder/dvf.py -i $f -o $HOME/project/DeepMicrobeFinder/result_other_2000/result_dvf/ -c 8
    # python $HOME/software/DeepVirFinder/dvf.py -i $f -o $HOME/project/DeepMicrobeFinder/data/virome/result_dvf/ -c 8
done