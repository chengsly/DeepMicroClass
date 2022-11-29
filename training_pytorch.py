import numpy as np
import constants
import sklearn
import optparse
import utils
import os, sys
from sklearn import preprocessing
from typing import Optional, Union
import argparse
import time

import torch
import torch.nn as nn
import torch.nn.functional as F
import pytorch_lightning as pl
from torch.utils.data import DataLoader, TensorDataset
from model.DeepMicroClass import DMF, DMFTransformer, LightningDMF, DMCLSTM
from pytorch_lightning import seed_everything, loggers
from sklearn.utils import class_weight
import pandas as pd
from model.SequenceData import SequenceDataset


####################################################################################################################################
#                                                      Set up                                                                      #
####################################################################################################################################
seed_everything(42)

parser = argparse.ArgumentParser()
parser.add_argument('--input', dest='input')
parser.add_argument('-l', '--log_prefix', dest='log_prefix', default='log')
parser.add_argument('--group', dest='group', default=None)

args = parser.parse_args()
log_prefix = args.log_prefix
group = args.group

if torch.cuda.is_available():
    device = torch.device("cuda")
    torch.backends.cudnn.benchmark = True
    print(torch.cuda.get_device_name(0))
else:
    device = torch.device("cpu")

TRAINING_DATASET_CSV = 'data/train.csv'
VAL_DATASET_CSV = 'data/val.csv'
ONEHOT_PATH = 'data/single_onehot'
BATCH_SIZE=128
LOG_DIR = f'data/pt_logs/{log_prefix}_{time.strftime("%Y-%m-%d-%H_%M", time.localtime())}'
CHECKPOINT_DIR = f'{LOG_DIR}/checkpoint'

# if not os.path.exists(CHECKPOINT_DIR):
#     os.makedirs(CHECKPOINT_DIR)

logger = loggers.WandbLogger(project='deepmicroclass', log_model=True, group=group)

############################################################################################################################
#                                             Training                                                                     #
############################################################################################################################


y = pd.read_table(TRAINING_DATASET_CSV, sep=',')
y = y['class']
weight = class_weight.compute_class_weight(class_weight='balanced', classes=np.unique(y), y=y)
training_weight = weight[y]
training_weight = torch.from_numpy(training_weight).float().to(device)

val_y = pd.read_table(VAL_DATASET_CSV, sep=',')
val_y = val_y['class']
val_weight = weight[val_y]
val_weight = torch.from_numpy(val_weight).float().to(device)

weight = torch.tensor(weight, dtype=torch.float).to(device)
print(f'Class weight: {weight}')

model = DMF()
# model = DMFTransformer()
# model = DMCLSTM()

logger.watch(model, log='all')

from pytorch_lightning.callbacks import ModelCheckpoint
trainer = pl.Trainer(
    accelerator='gpu',
    precision=16,
    max_epochs=100,
    default_root_dir=LOG_DIR,
    logger=logger,
    callbacks=[
        ModelCheckpoint(
            # dirpath=CHECKPOINT_DIR,
            filename='{epoch}-{step}-{val_f1:.3f}-{val_acc:.3f}',
            monitor='val_loss',
            # save_top_k=-1,
            mode='min',
            every_n_epochs=10, 
            save_last=True
            )
        ],
    # benchmark=True,
    )

train_dataset = SequenceDataset(TRAINING_DATASET_CSV, ONEHOT_PATH, preload=False)
val_dataset = SequenceDataset(VAL_DATASET_CSV, ONEHOT_PATH, preload=False)

train_dataloaders = DataLoader(
    train_dataset,
    batch_size=BATCH_SIZE,
    # shuffle=True,
    num_workers=16,
    collate_fn=SequenceDataset.custom_collate,
    sampler=torch.utils.data.sampler.WeightedRandomSampler(training_weight, BATCH_SIZE * 256)
)

val_dataloaders = DataLoader(
    val_dataset,
    batch_size=BATCH_SIZE,
    # shuffle=True,
    num_workers=16,
    collate_fn=SequenceDataset.custom_collate,
    sampler=torch.utils.data.sampler.WeightedRandomSampler(val_weight, BATCH_SIZE * 32)
)

trainer.fit(
    LightningDMF(model, weight=weight),
    train_dataloaders=train_dataloaders,
    val_dataloaders=val_dataloaders,
    # ckpt_path='data/pt_logs/checkpoint/epoch=1499-step=384000-val_f1=0.906-val_acc=0.907.ckpt'
)
