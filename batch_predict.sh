#!/bin/bash
#SBATCH --ntasks=1 --gres=gpu:1
#SBATCH --mem=64GB
#SBATCH --time=24:00:00

export PATH=$PATH:/home/rcf-40/siliangc/cmb/anaconda3/bin
source activate def
source /usr/usc/cuda/default/setup.sh
export LD_LIBRARY_PATH=/usr/usc/cuDNN/v7.6.5-cuda10.1/lib64:$LD_LIBRARY_PATH

time python predict_hybrid.py -i test.fasta -e one-hot -d models/one-hot-models/ -m single -l 500 
time python predict_hybrid.py -i test.fasta -e one-hot -d models/one-hot-models/ -m single -l 1000
time python predict_hybrid.py -i test.fasta -e one-hot -d models/one-hot-models/ -m single -l 2000
time python predict_hybrid.py -i test.fasta -e one-hot -d models/one-hot-models/ -m single -l 3000
time python predict_hybrid.py -i test.fasta -e one-hot -d models/one-hot-models/ -m single -l 5000
time python predict_hybrid.py -i test.fasta -e one-hot -d models/one-hot-models/ -m single
time python predict_hybrid.py -i test.fasta -e one-hot -d models/one-hot-models/ -m hybrid
