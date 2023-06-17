#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate pytorch

tensorboard --logdir_spec 500:500/lightning_logs,1000:1000/lightning_logs,2000:2000/lightning_logs,3000:3000/lightning_logs,5000:5000/lightning_logs