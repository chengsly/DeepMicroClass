#!/bin/sh

eval "$(conda shell.bash hook)"
# conda activate pytorch
conda activate pytorch

# bash batch_encode.sh data/training_data/ 500 one-hot
# bash batch_encode.sh data/training_data/ 1000 one-hot
# bash batch_encode.sh data/training_data/ 2000 one-hot
# bash batch_encode.sh data/training_data/ 3000 one-hot
# bash batch_encode.sh data/training_data/ 5000 one-hot

# python training.py -i data/local_training/encode_one-hot -e one-hot -t DNA -o data/new_model -l 500
# python training.py -i data/local_training/encode_one-hot -e one-hot -t DNA -o data/new_model -l 1000
# python training.py -i data/local_training/encode_one-hot -e one-hot -t DNA -o data/new_model -l 2000
# python training.py -i data/local_training/encode_one-hot -e one-hot -t DNA -o data/new_model -l 3000
# python training.py -i data/local_training/encode_one-hot -e one-hot -t DNA -o data/new_model -l 5000

python training_pytorch.py -i data/local_training/encode_one-hot -e one-hot -t DNA -o data/new_model -l 500
python training_pytorch.py -i data/local_training/encode_one-hot -e one-hot -t DNA -o data/new_model -l 1000
python training_pytorch.py -i data/local_training/encode_one-hot -e one-hot -t DNA -o data/new_model -l 2000
python training_pytorch.py -i data/local_training/encode_one-hot -e one-hot -t DNA -o data/new_model -l 3000
python training_pytorch.py -i data/local_training/encode_one-hot -e one-hot -t DNA -o data/new_model -l 5000

# python training.py -i data/training_data/encode_one-hot -e one-hot -t DNA -o data/new_model -l 500
# python training.py -i data/training_data/encode_one-hot -e one-hot -t DNA -o data/new_model -l 1000
# python training.py -i data/training_data/encode_one-hot -e one-hot -t DNA -o data/new_model -l 2000
# python training.py -i data/training_data/encode_one-hot -e one-hot -t DNA -o data/new_model -l 3000
# python training.py -i data/training_data/encode_one-hot -e one-hot -t DNA -o data/new_model -l 5000