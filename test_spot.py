import os, sys
from encoding_model import EncodingScheme
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
from scipy.special import softmax

import pytorch_lightning as pl
from model.DeepMicroClass import DeepMicroClass, LightningDMF
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
from model.SequenceData import SequenceDataset

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

model = LightningDMF.load_from_checkpoint('data/pt_logs/checkpoint/epoch=2999-step=768000-val_f1=0.906-val_acc=0.907.ckpt', model=DeepMicroClass(), map_location=device)
model.eval()
model.to(device)
model_cpu = LightningDMF.load_from_checkpoint('data/pt_logs/checkpoint/epoch=2999-step=768000-val_f1=0.906-val_acc=0.907.ckpt', model=DeepMicroClass())
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
name = 'Plants'
# in_fasta = 'data/SPOT/SPOT_5YrDuraVir_min_5000_newbler_toAmos_minimus2_id0.99_renamed.fa'
# in_fasta = '/home/shengwei/Viromes/04_read_mapping_per_sample/assemblies_per_sample/2SPOviralmetaG_derep.fasta'
in_fasta = f'/home/tianqi/dataset/IMGVR_v2/sources/{name}.fasta'
# in_fasta = '/home/tianqi/dataset/IMGVR_v4_high_confidence/sources/Human.fasta'



result_df = pd.DataFrame(columns=['id', 'Euk', 'EukVir', 'Plasmid', 'Prok', 'ProkVir'])

for record in tqdm(SeqIO.parse(in_fasta, 'fasta')):
    onehot = utils.seq2onehot(record.seq)
    contig_len = 2000
    # onehot = onehot[:onehot.shape[0]//5000*5000, :4].reshape(-1, 5000, 4)
    if onehot.shape[0] < contig_len:
        continue
    else:
        onehot = onehot[:onehot.shape[0]//contig_len*contig_len, :4].reshape(-1, contig_len, 4)
    onehot = onehot[:, None, :, :]
    onehot = torch.from_numpy(onehot).float()
    with torch.no_grad():
        try:
            pred = model(onehot.to(device))
            pred = pred.cpu().numpy()
        except:
            pred = model_cpu(onehot.to('cpu'))
            pred = pred.numpy()
    pred = softmax(pred, axis=1)
    pred = np.mean(pred, axis=0)
    result_df = pd.concat([result_df, pd.DataFrame([[record.id, pred[0], pred[1], pred[2], pred[3], pred[4]]], columns=['id', 'Euk', 'EukVir', 'Plasmid', 'Prok', 'ProkVir'])])

# out_path = 'data/SPOT/SPOT_5YrDuraVir_min_5000_newbler_toAmos_minimus2_id0.99_renamed_pred.csv'
out_path = f'data/virome/{name}.csv'
result_df.to_csv(out_path, index=False)

