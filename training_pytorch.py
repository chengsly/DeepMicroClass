import numpy as np
import constants
import sklearn
import optparse
import utils
import os, sys
from sklearn import preprocessing

import torch
import torch.nn as nn
import torch.nn.functional as F
import pytorch_lightning as pl
from torch.utils.data import DataLoader, TensorDataset
from DMF import DMF


####################################################################################################################################
#                                                      Set up                                                                      #
####################################################################################################################################

# POOL_FACTOR = 1
# dropout_cnn = 0.1
# dropout_pool = 0.1
# dropout_dense = 0.1
# learningrate = 0.001
# channel_num = 100    # 100-dim dense vector
# filters = [64, 128, 256]
# kernels = [6, 3, 2]
# ndense = 500
# epochs = 50

if torch.cuda.is_available():
    device = torch.device("cuda")
    torch.backends.cudnn.benchmark = True
    print(torch.cuda.get_device_name(0))
else:
    device = torch.device("cpu")
# torch.set_num_threads = 16
############################################################################################################################
#                                             Training                                                                     #
############################################################################################################################

from sklearn.utils import class_weight
import pandas as pd
# y = pd.read_table('data/seq_train.csv', sep=',')
# y = pd.read_table('data/local_training/train.csv', sep=',')
y = pd.read_table('data/train.csv', sep=',')
y = y['class']
# y = (y==0).astype(int)
weight = class_weight.compute_class_weight(class_weight='balanced', classes=np.unique(y), y=y)
# weight[2] *= 15
training_weight = weight[y]
training_weight = torch.from_numpy(training_weight).float().to(device)

# val_y = pd.read_table('data/local_training/test.csv', sep=',')
val_y = pd.read_table('data/test.csv', sep=',')
val_y = val_y['class']
val_weight = weight[val_y]
val_weight = torch.from_numpy(val_weight).float().to(device)

# weight = class_weight.compute_class_weight(class_weight='balanced', classes=np.unique(Y_tr_shuf.argmax(axis=1)), y=Y_tr_shuf.argmax(axis=1))
# weight = dict(enumerate(weight))
weight = torch.tensor(weight, dtype=torch.float).to(device)
print(weight)

batch_size=256
# pool_len1 = int((1000-500+1))

from DMF import LightningDMF, DMFTransformer

model = DMF()
# model = DMFTransformer()


from pytorch_lightning.callbacks import ModelCheckpoint, TQDMProgressBar
trainer = pl.Trainer(
    accelerator='gpu',
    precision=16,
    max_epochs=1500,
    # limit_train_batches=1024,
    default_root_dir=f'data/pt_logs/',
    callbacks=[
        ModelCheckpoint(dirpath=f'data/pt_logs/checkpoint_new/',
        filename='{epoch}-{step}-{val_f1:.3f}-{val_acc:.3f}',
        monitor='val_loss',
        save_top_k=-1, mode='min', every_n_epochs=10, save_last=True)
        ],
    # benchmark=True,
    )

# sampler = torch.utils.data.sampler.WeightedRandomSampler(weight, batch_size)

# train_dataloaders=DataLoader(TensorDataset(X_tr_shuf, Y_tr_shuf), batch_size=batch_size, shuffle=True, num_workers=8)
# val_dataloaders=DataLoader(TensorDataset(X_val, Y_val), batch_size=batch_size, num_workers=8)

# train_dataloaders=DataLoader(TensorDataset(X_tr_shuf, Y_tr_shuf), batch_size=batch_size, shuffle=True, num_workers=8)
from SequenceData import SequenceDataset
from torch.utils.data import default_collate
def custom_collate(batch):
    batch = list(filter(lambda x: x is not None, batch))
    # possible_contig_len = [500, 1000, 2000, 3000, 5000]
    # contig_len = possible_contig_len[np.random.randint(0, len(possible_contig_len))]
    contig_len = 5000
    batch = list((utils.sample_onehot(sample[0], contig_len)[None, :, :], sample[1]) for sample in batch)
    return default_collate(batch)

train_dataset = SequenceDataset('data/train.csv', 'data/single_onehot/', preload=False)
val_dataset = SequenceDataset('data/test.csv', 'data/single_onehot/', preload=False)

# train_dataset = SequenceDataset('data/local_training/train.csv', 'data/local_training/onehot/', preload=False)
# val_dataset = SequenceDataset('data/local_training/test.csv', 'data/local_training/onehot/', preload=False)

train_dataloaders = DataLoader(
    train_dataset,
    batch_size=batch_size,
    # shuffle=True,
    num_workers=16,
    collate_fn=custom_collate,
    sampler=torch.utils.data.sampler.WeightedRandomSampler(training_weight, batch_size * 256)
)

val_dataloaders = DataLoader(
    val_dataset,
    batch_size=batch_size,
    # shuffle=True,
    num_workers=16,
    collate_fn=custom_collate,
    sampler=torch.utils.data.sampler.WeightedRandomSampler(val_weight, batch_size * 32)
)

trainer.fit(
    LightningDMF(model, weight=weight),
    train_dataloaders=train_dataloaders,
    val_dataloaders=val_dataloaders,
    # ckpt_path='data/pt_logs/checkpoint/epoch=1499-step=384000-val_f1=0.906-val_acc=0.907.ckpt'
)
