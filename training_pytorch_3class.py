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
from torch.utils.data import DataLoader, TensorDataset, Dataset
from DMF import DMF_2class, DMF_3class

if torch.cuda.is_available():
    device = torch.device("cuda")
    torch.backends.cudnn.benchmark = True
    print(torch.cuda.get_device_name(0))
else:
    device = torch.device("cpu")

############################################################################################################################
#                                             Training                                                                     #
############################################################################################################################

class SequenceDataset(Dataset):
    def __init__(self, sequence_info, root_dir, preload=False):
        self.sequence_info = pd.read_table(sequence_info, sep=',',)
        self.sequence_info = self.sequence_info[(self.sequence_info.iloc[:, 1] == 2) | (self.sequence_info.iloc[:, 1] == 3) | (self.sequence_info.iloc[:, 1] == 4)]
        self.root_dir = root_dir
        self.preload = preload
        self.loading = True
        if preload:
            self.preloaded = []
            for i in range(len(self)):
                self.preloaded.append(self.__getitem__(i))
        self.loading = False

    def __len__(self):
        return len(self.sequence_info)

    def __getitem__(self, idx):
        id = self.sequence_info.iloc[idx, 0]
        class_ = self.sequence_info.iloc[idx, 1] - 2
        if self.preload and not self.loading:
            return self.preloaded[idx]
        else:
            sequence = np.load(os.path.join(self.root_dir, f'{id}.fasta.npy'))
            sequence = torch.from_numpy(sequence).float()
        return sequence, class_

class_for_train = 2

from sklearn.utils import class_weight
import pandas as pd
y = pd.read_table('data/train.csv', sep=',')
y = y[(y['class']==2)|(y['class']==3)|(y['class']==4)]
y = y['class'] - 2
weight = class_weight.compute_class_weight(class_weight='balanced', classes=np.unique(y), y=y)
weight[0] *= 10
training_weight = weight[y]
training_weight = torch.from_numpy(training_weight).float().to(device)

val_y = pd.read_table('data/test.csv', sep=',')
val_y = val_y[(val_y['class']==2)|(val_y['class']==3)|(val_y['class']==4)]
val_y = val_y['class'] - 2
val_weight = weight[val_y]
val_weight = torch.from_numpy(val_weight).float().to(device)

weight = torch.tensor(weight, dtype=torch.float).to(device)
print(weight)

batch_size=128
# pool_len1 = int((1000-500+1))

from DMF import LightningDMF

model = DMF_3class()

from pytorch_lightning.callbacks import ModelCheckpoint, TQDMProgressBar
trainer = pl.Trainer(
    accelerator='gpu',
    precision=16,
    max_epochs=1500,
    default_root_dir=f'data/pt_logs/',
    callbacks=[
        ModelCheckpoint(dirpath=f'data/pt_logs/checkpoint_prok/',
        filename='{epoch}-{step}-{val_loss:.2f}-{val_acc:.2f}',
        monitor='val_loss',
        save_top_k=-1, mode='min', every_n_epochs=10, save_last=True)
        ],
    # benchmark=True,
    )

from torch.utils.data import default_collate
def custom_collate(batch):
    batch = list(filter(lambda x: x is not None, batch))
    # possible_contig_len = [500, 1000, 2000, 3000, 5000]
    # contig_len = possible_contig_len[np.random.randint(0, len(possible_contig_len))]
    contig_len = 5000
    batch = list((utils.sample_onehot(sample[0], contig_len)[:, :], sample[1]) for sample in batch)
    return default_collate(batch)

train_dataset = SequenceDataset('data/train.csv', 'data/single_onehot/', preload=False)
val_dataset = SequenceDataset('data/test.csv', 'data/single_onehot/', preload=False)

train_dataloaders = DataLoader(
    train_dataset,
    batch_size=batch_size,
    shuffle=True,
    num_workers=8,
    collate_fn=custom_collate,
    # sampler=torch.utils.data.sampler.WeightedRandomSampler(training_weight, batch_size * 256)
)

val_dataloaders = DataLoader(
    val_dataset,
    batch_size=batch_size,
    shuffle=False,
    num_workers=8,
    collate_fn=custom_collate,
    # sampler=torch.utils.data.sampler.WeightedRandomSampler(val_weight, batch_size * 32)
)

trainer.fit(
    LightningDMF(model, weight=weight, num_classes=3),
    # train_dataloaders=DataLoader(TensorDataset(X_tr_shuf, Y_tr_shuf), batch_size=batch_size, num_workers=8, shuffle=True),
    train_dataloaders=train_dataloaders,
    # val_dataloaders=DataLoader(TensorDataset(X_val, Y_val), batch_size=batch_size, num_workers=8),
    val_dataloaders=val_dataloaders
)
