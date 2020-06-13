#!/bin/bash

#SBATCH --ntasks=1 --gres=gpu:2
#SBATCH --mem=181GB
#SBATCH --time=24:00:00

export PATH=$PATH:/home/rcf-40/siliangc/cmb/anaconda3/bin
source activate dvf
source /usr/usc/cuda/default/setup.sh
export LD_LIBRARY_PATH=/usr/usc/cuDNN/v7.6.5-cuda10.1/lib64:$LD_LIBRARY_PATH
python training_lstm.py -i $1 -e $2 -l $3 -t $4 -o $5 
