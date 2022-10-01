from torch.utils.data import Dataset, DataLoader
import torch
import pandas as pd
import os
import numpy as np

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