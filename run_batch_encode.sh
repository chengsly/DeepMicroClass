#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate tf2

bash batch_encode.sh data/training_data/ 500 one-hot
bash batch_encode.sh data/training_data/ 1000 one-hot
bash batch_encode.sh data/training_data/ 2000 one-hot
bash batch_encode.sh data/training_data/ 3000 one-hot
bash batch_encode.sh data/training_data/ 5000 one-hot