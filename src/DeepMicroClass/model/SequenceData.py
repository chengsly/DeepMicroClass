from torch.utils.data import Dataset, DataLoader
import torch
import pandas as pd
import os
import numpy as np
from typing import Optional, Union
from torch.utils.data import default_collate

import utils

class SequenceDataset(Dataset):
    def __init__(self, sequence_info, root_dir, preload=False):
        self.sequence_info = pd.read_table(sequence_info, sep=',',)
        self.c = np.zeros(self.sequence_info.shape[0])
        self.c[self.sequence_info.iloc[:, 1]==0] = 1
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
        class_ = self.sequence_info.iloc[idx, 1]
        # class_ = self.c[idx]
        if self.preload and not self.loading:
            return self.preloaded[idx]
        else:
            sequence = np.load(os.path.join(self.root_dir, f'{id}.fasta.npy'))
            sequence = torch.from_numpy(sequence).float()
        return sequence, class_

    @staticmethod
    def custom_collate(batch, contig_len:Optional[Union[list, int]]=None):
        batch = list(filter(lambda x: x is not None, batch))
        # possible_contig_len = [500, 1000, 2000, 3000, 5000]
        # contig_len = possible_contig_len[np.random.randint(0, len(possible_contig_len))]
        contig_len = 5000
        batch = list((utils.sample_onehot(sample[0], contig_len)[None, :, :], sample[1]) for sample in batch)
        return default_collate(batch)

class SequenceDataset_tfidf(Dataset):
    def __init__(self, sequence_info, root_dir, tfidf_dir, preload=False):
        self.sequence_info = pd.read_table(sequence_info, sep=',',)
        self.c = np.zeros(self.sequence_info.shape[0])
        self.c[self.sequence_info.iloc[:, 1]==0] = 1
        self.root_dir = root_dir
        self.tfidf_dir = tfidf_dir
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
        # class_ = self.sequence_info.iloc[idx, 1]
        class_ = self.c[idx]
        if self.preload and not self.loading:
            return self.preloaded[idx]
        else:
            sequence = np.load(os.path.join(self.root_dir, f'{id}.fasta.npy'))
            sequence = torch.from_numpy(sequence).float()
            tfidf = torch.from_numpy(np.load(os.path.join(self.tfidf_dir, f'{id}.fasta.npy'))).float()
        return sequence, tfidf, class_