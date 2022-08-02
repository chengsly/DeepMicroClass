from Bio import SeqIO
import pandas as pd
import os
from ete3 import NCBITaxa
import numpy as np
import math
from random import sample

PROK_PATH = "data/filtered_fasta/prokaryote.fa"
EUK_PATH = "data/filtered_fasta/eukaryote.fa"
PLASMID_PATH = "data/filtered_fasta/plasmid.fa"
PROKVIR_PATH = "data/filtered_fasta/prok_vir.fa"
EUKVIR_PATH = "data/filtered_fasta/euk_vir.fa"

PATHS = [PROK_PATH, EUK_PATH, PLASMID_PATH, PROKVIR_PATH, EUKVIR_PATH]

TARGET_PATH = "data/sampled_sequence"
os.makedirs(TARGET_PATH, exist_ok=True)

TOTAL_SEQ_NUM = 1000

PROK_EUK_RATIO = np.array([
    [9, 1],
    [7, 3],
    [5, 5],
    [3, 7],
    [1, 9],
], dtype=float)
PROK_EUK_RATIO /= PROK_EUK_RATIO.sum(axis=1, keepdims=True)

PROK_PROKVIR_PLASMID_RATIO = np.array([
    [5, 1, 1],
    [4, 1, 1],
    [3, 1, 1],
    [2, 1, 1],
], dtype=float)
PROK_PROKVIR_PLASMID_RATIO /= PROK_PROKVIR_PLASMID_RATIO.sum(axis=1, keepdims=True)

EUK_EUKVIR_RATIO = np.array([
    [5, 1],
    [4, 1],
    [3, 1],
    [2, 1],
], dtype=float)
EUK_EUKVIR_RATIO /= EUK_EUKVIR_RATIO.sum(axis=1, keepdims=True)

# for i in range(len(PROK_EUK_RATIO)):
#     for j in range(len(PROK_PROKVIR_PLASMID_RATIO)):

#         # prok_total = math.floor(PROK_EUK_RATIO[i, 0] * TOTAL_SEQ_NUM)
#         # euk_total = math.floor(PROK_EUK_RATIO[i, 1] * TOTAL_SEQ_NUM)

ratio = np.concatenate([np.concatenate([PROK_EUK_RATIO[i, 0] * PROK_PROKVIR_PLASMID_RATIO, PROK_EUK_RATIO[i, 1] * EUK_EUKVIR_RATIO], axis=1) for i in range(PROK_EUK_RATIO.shape[0])], axis=0)
seq_num = np.round(ratio * TOTAL_SEQ_NUM).astype(int)
print(seq_num)

records = [list(SeqIO.parse(path, 'fasta')) for path in PATHS]

for i in range(seq_num.shape[0]):
    sampled_records = [sample(records[j], seq_num[i, j]) for j in range(seq_num.shape[1])]
    output_path = os.path.join(TARGET_PATH, f'{i:02d}_PROK{len(sampled_records[0])}_PROKVIR{len(sampled_records[1])}_PLSMD{len(sampled_records[2])}_EUK{len(sampled_records[3])}_EUKVIR{len(sampled_records[4])}.fa')
    SeqIO.write(sum(sampled_records, []), output_path, 'fasta')