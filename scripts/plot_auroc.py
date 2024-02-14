import os, sys
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np

import torch
import pickle
from sklearn.metrics import roc_auc_score, f1_score, accuracy_score
from torchmetrics.functional import auroc, roc, auc
import matplotlib.pyplot as plt
from tqdm import tqdm

contig_lengths = [500, 1000, 2000, 3000, 5000, 10000, 50000, 100000]
# contig_lengths = [500, 1000, 2000, 3000, 5000, 10000, 50000]
labels = ['Euk', 'EukVir', 'Plasmid', 'Prok', 'ProkVir']

fig, axes = plt.subplots(4, 2, figsize=(8, 10), sharex=True, sharey=True)

for i, contig_len in enumerate(contig_lengths):
    fn = f'results/contig_len_{contig_len}.pkl'
    with open(fn, 'rb') as f:
        data = pickle.load(f)
        all_y = data['all_y']
        all_y_hat = data['all_y_hat']
    # print(auroc(torch.tensor(all_y_hat), torch.tensor(all_y), num_classes=5, average='None'))
    fpr, tpr, _ = roc(torch.tensor(all_y_hat), torch.tensor(all_y), num_classes=5)
    for j in range(5):
        axes[i//2, i%2].plot(fpr[j].numpy(), tpr[j].numpy(), label=f'{labels[j]}, AUC={auc(fpr[j], tpr[j]):.3f}', linewidth=0.5)
    axes[i//2, i%2].legend(prop={'size': 9})
    axes[i//2, i%2].set_title(f'Contig length: {contig_len} bps', fontsize=8)
plt.yticks(np.arange(0, 1.01, 0.2))
# plt.xlabel('False Positive Rate')
fig.supxlabel('False Positive Rate', fontsize=16)
fig.supylabel('True Positive Rate', fontsize=16)
fig.tight_layout(pad=0.8)
# for i, contig_len in enumerate(contig_lengths):
#     fn = f'results/contig_len_{contig_len}.pkl'
#     with open(fn, 'rb') as f:
#         data = pickle.load(f)
#         all_y = data['all_y']
#         all_y_hat = data['all_y_hat']

#     print(auroc(torch.tensor(all_y_hat), torch.tensor(all_y), num_classes=5, average='None'))
#     fpr, tpr, _ = roc(torch.tensor(all_y_hat), torch.tensor(all_y), num_classes=5)
#     for j in range(5):
#         axes[i].plot(fpr[j].numpy(), tpr[j].numpy(), label=f'{labels[j]}')
#     axes[i].legend(prop={'size': 3})
plt.savefig('auroc.pdf', bbox_inches='tight')