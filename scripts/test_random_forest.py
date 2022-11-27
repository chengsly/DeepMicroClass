import numpy as np
import constants
import sklearn
import optparse
import utils
import os, sys
from sklearn import preprocessing
from typing import Union
import argparse
import time

import torch
import torch.nn as nn
import torch.nn.functional as F
import pytorch_lightning as pl
from torch.utils.data import DataLoader, TensorDataset
from pytorch_lightning import seed_everything, loggers
from sklearn.ensemble import RandomForestClassifier
from sklearn.utils import class_weight
import pandas as pd
from SequenceData import SequenceDataset

import wandb

seed_everything(42)

TRAINING_DATASET_CSV = 'data/train.csv'
VAL_DATASET_CSV = 'data/val.csv'
ONEHOT_PATH = 'data/single_onehot'
BATCH_SIZE=128

wandb.init(project='deepmicroclass', entity='deepmicroclass', name='random_forest')

y = pd.read_table(TRAINING_DATASET_CSV, sep=',')
y = y['class']
weight = class_weight.compute_class_weight(class_weight='balanced', classes=np.unique(y), y=y)
training_weight = weight[y]

val_y = pd.read_table(VAL_DATASET_CSV, sep=',')
val_y = val_y['class']
val_weight = weight[val_y]

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

model = RandomForestClassifier(n_estimators=100, max_depth=10, class_weight='balanced')