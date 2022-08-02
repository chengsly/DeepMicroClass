import os, sys
from encoding_model import EncodingScheme
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import constants
import multiprocessing

import pytorch_lightning as pl
from DMF import DMF, LightningDMF
import torch
import torch.nn as nn
import torch.nn.functional as F
import pandas as pd
from functools import partial
from torch.utils.data import DataLoader, Subset
import pickle
from sklearn.metrics import roc_auc_score, f1_score, accuracy_score
from torchmetrics.functional import auroc, roc
import matplotlib.pyplot as plt
from tqdm import tqdm
import copy

import utils
from SequenceData import SequenceDataset

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

model = LightningDMF.load_from_checkpoint('data/pt_logs/checkpoint/epoch=2999-step=768000-val_f1=0.906-val_acc=0.907.ckpt', model=DMF(), map_location=device)
model.eval()
model.to(device)
model_cpu = LightningDMF.load_from_checkpoint('data/pt_logs/checkpoint/epoch=2999-step=768000-val_f1=0.906-val_acc=0.907.ckpt', model=DMF())
model_cpu.eval()
batch_size = 32

from torch.utils.data import default_collate
def custom_collate(batch, contig_len=5000):
    batch = list(filter(lambda x: x is not None, batch))
    batch = list(filter(lambda x: x[0].shape[0] >= contig_len, batch))
    # possible_contig_len = [500, 1000, 2000, 3000, 5000]
    # contig_len = possible_contig_len[np.random.randint(0, len(possible_contig_len))]
    # contig_len = contig_len
    batch = list((utils.sample_onehot(sample[0], contig_len)[None, :, :], sample[1]) for sample in batch)
    if len(batch) == 0:
        return None
    return default_collate(batch)

test_dataset = SequenceDataset('data/final_test.csv', 'data/single_onehot/', preload=False)

# indices = torch.randperm(len(test_dataset))[:3000]
# test_test_dataset = Subset(test_dataset, indices)

# contig_lengths = [500, 1000, 2000, 3000, 5000]
contig_lengths = [500, 1000, 2000, 3000, 5000, 10000, 50000, 100000]

torch.set_num_threads(16)

with torch.no_grad():
    for contig_len in contig_lengths:
        test_dataloaders = DataLoader(
            test_dataset,
            batch_size=batch_size,
            shuffle=False,
            num_workers=8,
            collate_fn=partial(custom_collate, contig_len=contig_len),
        )
        # trainer = pl.Trainer(gpus=1, precision=16)
        # metrics = trainer.test(model, test_dataloaders)
        # print(metrics)
        # with open(f'results/contig_len_{contig_len}.pkl', 'wb') as f:
        #     pickle.dump(metrics, f)
        fn = f'results/contig_len_{contig_len}.pkl'
        if os.path.exists(fn):
            with open(fn, 'rb') as f:
                data = pickle.load(f)
                all_y = data['all_y']
                all_y_hat = data['all_y_hat']
        else:
            all_y_hat = []
            all_y = []
            for step, batch in enumerate(tqdm(test_dataloaders)):
                # print(step, end='\r')
                if batch is None:
                    continue
                x, y = batch
                x = x.to(device)
                y = y.to(device)
                y_hat = model(x)
                all_y_hat.append(y_hat.cpu().numpy())
                all_y.append(y.cpu().numpy())
            all_y_hat = np.concatenate(all_y_hat)
            all_y = np.concatenate(all_y)
            with open(fn, 'wb') as f:
                pickle.dump({'all_y': all_y, 'all_y_hat': all_y_hat}, f)
        # print(roc_auc_score(all_y, all_y_hat))
        print(auroc(torch.tensor(all_y_hat), torch.tensor(all_y), num_classes=5, average='None'))
        print((all_y == 0).sum(), (all_y == 1).sum(), (all_y == 2).sum(), (all_y == 3).sum(), (all_y == 4).sum())
    
    # fn = f'results/whole.pkl'
    # if os.path.exists(fn):
    #     with open(fn, 'rb') as f:
    #         data = pickle.load(f)
    #         all_y = data['all_y']
    #         all_y_hat = data['all_y_hat']
    #         print(auroc(torch.tensor(all_y_hat), torch.tensor(all_y), num_classes=5, average='None'))
    # else:
    #     all_y_hat = []
    #     all_y = []
    #     for x, y in tqdm(test_dataset):
    #         x = x[None, None, :, :]
    #         try:
    #             x.to(device)
    #             y_hat = model(x)
    #             all_y_hat.append(y_hat.cpu().numpy())
    #             all_y.append(y)
    #         except:
    #             x.to('cpu')
    #             y_hat = model_cpu(x)
    #             all_y_hat.append(y_hat.numpy())
    #             all_y.append(y)
    #     all_y_hat = np.concatenate(all_y_hat)
    #     all_y = np.concatenate(all_y)
    #     with open(fn, 'wb') as f:
    #         pickle.dump({'all_y': all_y, 'all_y_hat': all_y_hat}, f)

