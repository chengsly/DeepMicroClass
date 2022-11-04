#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate ppr-meta

export LD_LIBRARY_PATH=/home/tianqi/software/MATLABRuntime/v94/runtime/glnxa64:/home/tianqi/software/MATLABRuntime/v94/bin/glnxa64:/home/tianqi/software/MATLABRuntime/v94/sys/os/glnxa64:/home/tianqi/software/MATLABRuntime/v94/extern/bin/glnxa64

for f in /home/tianqi/project/DeepMicrobeFinder/data/sampled_sequence/*.fa
do
    filename=$(basename -- "$f")
    extension="${filename##*.}"
    filename="${filename%.*}"
    echo "Predicting $f"
    /home/tianqi/software/PPR-Meta/PPR_Meta $f result/${filename}.csv
done
