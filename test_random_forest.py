import numpy as np
import sklearn
import optparse
import os, sys
# import utils
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

from model.SequenceData import SequenceDataset
from model.DeepMicroClass import KMerTransformer

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
    sampler=torch.utils.data.sampler.WeightedRandomSampler(val_weight, BATCH_SIZE * 256)
)

# transformer = KMerTransformer(k=3, rearrange=True)

train_X = []
train_y = []
for batch in train_dataloaders:
    X = batch[0]
    # X = transformer(X)
    # X = torch.argmax(X, dim=1).long()
    train_X.append(X.numpy())
    train_y.append(batch[1].numpy())
train_X = np.concatenate(train_X, axis=0)
train_y = np.concatenate(train_y, axis=0)
train_X = train_X.reshape(train_X.shape[0], -1)

val_X = []
val_y = []
for batch in val_dataloaders:
    X = batch[0]
    # X = transformer(X)
    # X = torch.argmax(X, dim=1).long()
    val_X.append(X.numpy())
    val_y.append(batch[1].numpy())
val_X = np.concatenate(val_X, axis=0)
val_y = np.concatenate(val_y, axis=0)
val_X = val_X.reshape(val_X.shape[0], -1)

model = RandomForestClassifier(n_estimators=100, max_depth=10, class_weight='balanced', n_jobs=16)

model.fit(train_X, train_y)
pred = model.predict(val_X)
"""Calculate validation accuracy and F1 score, and log them to wandb"""
val_acc = (pred == val_y).sum() / len(val_y)
val_f1 = sklearn.metrics.f1_score(val_y, pred, average='weighted')
wandb.log({'val_acc': val_acc, 'val_f1': val_f1})