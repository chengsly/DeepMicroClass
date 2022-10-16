 1/1: import pandas as pd
 1/2:
DB_METADATA_PATH = "data/IMGVR/IMGVR_all_Host_information.tsv"
metadata = pd.read_table(DB_METADATA_PATH)
 1/3: metadata.columns
 1/4: metadata['Host taxonomy prediction'][0]
 1/5:
host = set()
for _, row in metadata.iterrows():
    host.add(row['Host taxonomy prediction'].split(';')[0])
 1/6: host
 1/7:
FASTA_PATH = "data/IMGVR/IMGVR_all_nucleotides.fna"
from Bio import SeqIO
a = next(SeqIO.parse(FASTA_PATH, "fasta"))
 1/8: a
 1/9:
DB_METADATA_PATH = "data/IMGVR/IMGVR_all_Sequence_information.tsv"
metadata = pd.read_table(DB_METADATA_PATH)
1/10: metadata.columns
1/11:
host = set()
for _, row in metadata.iterrows():
    host.add(row['Host taxonomy prediction'].split(';')[0])
1/12: metadata.columns[17]
1/13: a = metadata['Host taxonomy prediction']
1/14:
metadata.dropna(subset=["Host taxonomy prediction"], inplace=True)
metadata.reset_index(drop=True, inplace=True)
1/15:
host = set()
for _, row in metadata.iterrows():
    host.add(row['Host taxonomy prediction'].split(';')[0])
1/16: host
1/17: a.id
1/18: a = next(SeqIO.parse(FASTA_PATH, "fasta"))
1/19: a.id.split('|')
 2/1:
import pandas as pd
dist = pd.read_table('data/dist/all_prokaryote.tsv')
 3/1:
import pandas as pd
dist = pd.read_table('data/dist/refseq_eukaryote.tsv')
 3/2: dist = pd.read_table('data/dist/refseq_eukaryote.tsv', index_col=0)
 3/3:
from Bio import SeqIO
r = next(SeqIO.parse('data/refseq_eukaryote.fa'))
 3/4:
from Bio import SeqIO
r = next(SeqIO.parse('data/refseq_eukaryote.fa', 'fasta'))
 3/5: r
 3/6: dist[r.id, r.id]
 3/7: dist[r.id][r.id]
 3/8:
import numpy as np
a = np.load('/home/tianqi/project/siliangc/DeepEukFinder/allSeqs/ProkVal#ProkVal#2k_seq10383_codefw.npy')
 3/9: dist.min()
3/10: dist[dist>0].min()
3/11: dist[['FR796397.1', 'FR796398.1']]
3/12: dist[['FR796397.1', 'FR796398.1']][['FR796397.1', 'FR796398.1']]
3/13: dist[['FR796397.1', 'FR796398.1']][['FR796397.1', 'FR796398.1'],:]
3/14: dist.loc[['FR796397.1', 'FR796398.1'], ['FR796397.1', 'FR796398.1']]
3/15: b = np.arange(10)
3/16: c = [i for i in b]
3/17: dist.loc[['FR796397.1', 'FR796398.1'], ['FR796397.1']]
3/18: dist.loc[['FR796397.1', 'FR796398.1'], ['FR796397.1']]<0.2
3/19: c = dist.loc[['FR796397.1', 'FR796398.1'], ['FR796397.1']]
3/20: c.where(c<0.01)
3/21: c.index(c<0.01)
3/22: c.index[c<0.01]
3/23: c<0.01
3/24: (c<0.01).tolist()
3/25: c = dist.loc[['FR796397.1'], ['FR796397.1', 'FR796398.1']]
3/26: c
3/27: c.index[c<0.01]
3/28: _, row = next(c.iterrows())
3/29: row
3/30: row<0.1
3/31: row.index[row<0.1]
3/32: c.index
3/33: c = dist.loc[['FR796397.1', 'FR796398.1'], ['FR796397.1', 'FR796398.1']]
3/34: c
3/35: c.index
3/36: c.index[c<0.01]
3/37: c[c<0.01].sum(axis=1)
3/38: c[c<0.01].select_dtypes(include=['bool']).sum(axis=1)
3/39: c[c<0.01]
3/40: (c<0.01).select_dtypes(include=['bool']).sum(axis=1)
 4/1:
import numpy as np
a = np.load('data/training_data/encode_one-hot/EukVirusVal#EukVirusVal#0.5k_seq17154_codebw.npy')
 4/2: np.savez_compressed('data/training_data/encode_one-hot/EukVirusVal#EukVirusVal#0.5k_seq17154_codebw.npy', a=a)
 4/3: np.save('data/training_data/encode_one-hot/EukVirusVal#EukVirusVal#0.5k_seq17154_codebw.npy', a)
 4/4: import gzip
 4/5:
b = gzip.open('data/training_data/encode_one-hot/EukVirusVal#EukVirusVal#0.5k_seq17154_codebw.npy.gz', mode='wt')
np.save(b, a)
 4/6:
b = gzip.open('data/training_data/encode_one-hot/EukVirusVal#EukVirusVal#0.5k_seq17154_codebw.npy.gz', mode='w')
np.save(b, a)
 4/7:
import pgzip
b = pgzip.open('data/training_data/encode_one-hot/EukVirusVal#EukVirusVal#0.5k_seq17154_codebw.npy.gz', mode='w', thread=16)
np.save(b, a)
 4/8:
b = pgzip.open('data/training_data/encode_one-hot/EukVirusVal#EukVirusVal#0.5k_seq17154_codefw.npy.gz', thread=16)
c = np.load(b)
 4/9: np.array_equal(a, c)
4/10: a.dtype
4/11: c.dtype
4/12: np.array_equal(a, c[::-1, ::-1])
4/13: np.array_equal(a, c[,::-1, ::-1])
4/14: np.array_equal(a, c[:,::-1, ::-1])
4/15: b.close()
 5/1:
import pandas as pd
df = pd.read_table('data/org_Archaea/assembly_summary.txt')
 5/2: df = pd.read_table('data/org_Archaea/assembly_summary.txt', column=None)
 5/3: df = pd.read_table('data/org_Archaea/assembly_summary.txt', header=None)
 5/4: df[:, 14]
 5/5: a = df[14].to_datetime()
 5/6: a = pd.to_datetime(df[14])
 5/7: import datetime
 5/8: datetime(2020, 1, 1) < a
 5/9: b = datetime.datetime(2020, 1, 1) < a
5/10:
import numpy as np
b = np.datetime64('2020-01-01') < a
5/11: b[348]
5/12: df = pd.read_table('data/org_Archaea/assembly_summary.txt')
5/13: df = pd.read_table('data/org_Archaea/assembly_summary.txt', header=None)
5/14: df[6].value_counts()
5/15: a = df[6].value_counts()
5/16: df = pd.read_table('data/org_Bacteria/assembly_summary.txt', header=None)
5/17: a = df[6].value_counts()
5/18: df = pd.read_table('data/pre20_bacteria/assembly_summary.txt', header=None)
5/19: a = df[6].value_counts()
5/20:
df = pd.read_table('data/post20_bacteria/assembly_summary.txt', header=None)
a = df[6].value_counts()
5/21:
b = np.arange(10)
b[:20]
5/22: df = pd.read_table('data/plsdb/plsdb.tsv')
5/23: dates = df['SubmissionDate_ASSEMBLY']
5/24: dates = pd.to_datetime(dates)
5/25: b = np.datetime64('2020-01-01')>dates
5/26: b[2]
5/27: pre20 = df.loc[b, :]
5/28: b = dates < np.datetime64('2020-01-01')
 6/1:
import pandas as pd
df = pd.read_tabel('DeepEukFinder_resultst/eukaryote.fa_pred_one-hot_hybrid.txt')
 6/2: df = pd.read_table('DeepEukFinder_resultst/eukaryote.fa_pred_one-hot_hybrid.txt')
 6/3: result = df.iloc[:, 1:]
 6/4: result_max = result.argmax(axis=1)
 6/5:
import numpy as np
result = result.to_numpy()
 6/6: result_max = result.argmax(axis=1)
 6/7: (result_max==0).sum()/result_max.shape
 6/8:
from Bio import SeqIO
a = next(SeqIO.parse('data/euk_fasta/GCA_003072545.2_ASM307254v2_genomic.fna', 'fasta'))
 6/9: b = str(a)
6/10: b = str(a.seq)
6/11: b = str(a.seq.seq)
 7/1:
import pandas as pd
result = pd.read_table('DeepMicrobeFinder_results/00_PROK643_PROKVIR129_PLSMD129_EUK83_EUKVIR17.fa_pred_one-hot_hybrid.txt')
 7/2: result_np = result.iloc[:, 1:].to_numpy()
 7/3: max_idx = result_np.argmax(axis=1)
 7/4: (max_idx[:643]==3).sum()/643
 7/5: (max_idx[643:643+129]==4).sum()/129
 7/6: (max_idx[643+129:643+129+129]==2).sum()/129
 7/7: (max_idx[643+129+129:643+129+129+83]==2).sum()/83
 7/8: (max_idx[643+129+129:643+129+129+83]==0).sum()/83
 7/9: (max_idx[643+129+129+83]==0).sum()/17
7/10: (max_idx[643+129+129+83]==1).sum()/17
7/11: (max_idx[643+129+129+83:]==1).sum()/17
7/12: target = np.array([3]*643+[4]*129+[2]*129+[0]*83+[1]*17)
7/13:
import numpy as np
target = np.array([3]*643+[4]*129+[2]*129+[0]*83+[1]*17)
7/14: (max_idx == target).sum()/1001
7/15: target = np.array([3]*643+[0]*129+[2]*129+[4]*83+[1]*17)
7/16: (max_idx == target).sum()/1001
7/17: (max_idx[643+129+129:643+129+129+83]==4).sum()/83
7/18: (max_idx[643:643+129]==0).sum()/129
7/19: a = pd.read_table('data/result_tiara/00_PROK643_EUK129_PLSMD129_PROKVIR83_EUKVIR17.txt')
7/20:
from Bio import SeqIO
seqs = list(SeqIO.parse('DeepMicrobeFinder_results/19_PROK50_EUK25_PLSMD25_PROKVIR600_EUKVIR300.fa_pred_one-hot_hybrid.txt', 'fasta'))
7/21:
from Bio import SeqIO
seqs = list(SeqIO.parse('data/sampled_sequence/00_PROK643_EUK129_PLSMD129_PROKVIR83_EUKVIR17.fa', 'fasta'))
7/22: seqs[0]
7/23: ids = [seq.id for seq in seqs]
7/24:
def construct_target(nums):
    order = [[3], [0], [2], [4], [1]]
    return np.array(sum([order[i]*nums[i] for i in range(5)], []))
7/25:
import re
f = '00_PROK643_EUK129_PLSMD129_PROKVIR83_EUKVIR17.fa'
nums = re.findall(r'\d+', f)[1:]
nums = [int(n) for n in nums]
target = construct_target(nums)
7/26: df = pd.DataFrame(target, index=ids)
7/27: a = pd.read_table('data/result_tiara/00_PROK643_PROKVIR129_PLSMD129_EUK83_EUKVIR17.txt')
7/28: a = pd.read_table('data/result_tiara/00_PROK643_PROKVIR129_PLSMD129_EUK83_EUKVIR17.txt', index_col=0)
7/29: a.index
7/30: a.index[0]
7/31: a.index[0].split()
 8/1: a = 1
 9/1:
import pandas as pd
import numpy as np
10/1:
import pandas as pd
import numpy as np
10/2:
import numpy as np
import pandas as pd
import os
import re
from sklearn.metrics import f1_score
10/3:
RESULT_DIR = 'data/result_tiara'
results_fn = os.listdir(RESULT_DIR)
results_fn = [f for f in results_fn if not f.startswith('log')]
results_fn.sort()
10/4:
result = pd.read_table(os.path.join(RESULT_DIR, f), index_col=0)
index = list(result.index)
index = [i.split('')[0] for i in index]
result.index = index
10/5:
result = pd.read_table(os.path.join(RESULT_DIR, results_fn[0]), index_col=0)
    
index = list(result.index)
index = [i.split('')[0] for i in index]
result.index = index
10/6:
index = list(result.index)
index = [i.split()[0] for i in index]
result.index = index
10/7: target = pd.read_pickle('data/sampled_sequence_target/00_PROK643_PROKVIR129_PLSMD129_EUK83_EUKVIR17.fa.pkl')
10/8: actual_target = [target[i] for i in index]
10/9: actual_target = [target.loc[i] for i in index]
10/10: a = pd.read_pickle('data/sampled_sequence_target/03_PROK450_PROKVIR225_PLSMD225_EUK67_EUKVIR33.fa.pkl')
10/11: a = pd.read_pickle('data/sampled_sequence_target/02_PROK540_PROKVIR180_PLSMD180_EUK75_EUKVIR25.fa.pkl')
10/12: a = pd.read_pickle('data/sampled_sequence_target/02_PROK540_PROKVIR180_PLSMD180_EUK75_EUKVIR25.fa.pkl')
10/13: a = pd.read_pickle('data/sampled_sequence_target/02_PROK540_PROKVIR180_PLSMD180_EUK75_EUKVIR25.fa.pkl')
10/14: a = pd.read_pickle('data/sampled_sequence_target/03_PROK450_PROKVIR225_PLSMD225_EUK67_EUKVIR33.fa.pkl')
10/15: a = pd.read_pickle('data/sampled_sequence_target/03_PROK450_PROKVIR225_PLSMD225_EUK67_EUKVIR33.fa.pkl')
11/1:
import pandas as pd
import numpy as np
a = pd.read_pickle('data/sampled_sequence_target/03_PROK450_PROKVIR225_PLSMD225_EUK67_EUKVIR33.fa.pkl')
12/1:
from Bio import SeqIO
prok = list(SeqIO.parse('data/org_fasta/vhdb_prok_vir.fa'), 'fasta')
euk = list(SeqIO.parse('data/org_fasta/vhdb_euk_vir.fa'), 'fasta')
12/2:
from Bio import SeqIO
prok = list(SeqIO.parse('data/org_fasta/vhdb_prok_vir.fa', 'fasta'))
euk = list(SeqIO.parse('data/org_fasta/vhdb_euk_vir.fa', 'fasta'))
12/3:
pid = [seq.id for seq in prok]
eid = [seq.id for seq in euk]
12/4: set(pid) & set(eid)
12/5: len(set(pid) | set(eid))
12/6:
import pandas as pd
df = pd.read_pickle('data/sampled_sequence_target/03_PROK450_PROKVIR225_PLSMD225_EUK67_EUKVIR33.fa.pkl')
12/7: df.iloc[117, :]
12/8:
del prok
del euk
12/9:
bac = [seq.id for seq in SeqIO.parse('data/org_fasta/post20_bacteria.fa', 'fasta')]
pls = [seq.id for seq in SeqIO.parse('data/plasmid/post20_plasmid.fa', 'fasta')]
12/10: set(bac) & set(pls)
12/11: len(set(bac) & set(pls))
13/1: import pandas as pd
13/2: df = pd.read_table('data/result_whokaryote/00_PROK643_PROKVIR129_PLSMD129_EUK83_EUKVIR17/whokaryote_predictions_T.tsv')
13/3: df = pd.read_table('data/result_whokaryote/00_PROK643_PROKVIR129_PLSMD129_EUK83_EUKVIR17/whokaryote_predictions_T.tsv', index_col=0)
13/4: target = pd.read_table('data/sampled_sequence_target/00_PROK643_PROKVIR129_PLSMD129_EUK83_EUKVIR17.fa.pkl')
13/5: target = pd.read_pickle('data/sampled_sequence_target/00_PROK643_PROKVIR129_PLSMD129_EUK83_EUKVIR17.fa.pkl')
14/1: %run scripts/parse_tiara_output.py
14/2: %run scripts/parse_tiara_output.py
14/3: %run scripts/parse_tiara_output.py
14/4: %run scripts/parse_tiara_output.py
14/5: %run scripts/parse_tiara_output.py
14/6: %run scripts/parse_tiara_output.py
14/7: %run scripts/parse_tiara_output.py
14/8: %run scripts/parse_tiara_output.py
14/9: %run scripts/parse_tiara_output.py
14/10: %run scripts/parse_tiara_output.py
14/11: %run scripts/parse_tiara_output.py
14/12: %run scripts/parse_tiara_output.py
14/13: %run scripts/parse_tiara_output.py
14/14: %run scripts/parse_tiara_output.py
15/1: %run training_pytorch.py
15/2: %run training_pytorch.py
15/3: %run training_pytorch.py -i data/local_training/encode_one-hot -e one-hot -t DNA -o data/new_model -l 500
15/4: %run training_pytorch.py -i data/local_training/encode_one-hot -e one-hot -t DNA -o data/new_model -l 500
15/5: %run training_pytorch.py -i data/local_training/encode_one-hot -e one-hot -t DNA -o data/new_model -l 500
15/6:
%load_ext autoreload
%autoreload 2
15/7: %run training_pytorch.py -i data/local_training/encode_one-hot -e one-hot -t DNA -o data/new_model -l 500
16/1:
%load_ext autoreload
%autoreload 2
16/2: %run training_pytorch.py -i data/local_training/encode_one-hot -e one-hot -t DNA -o data/new_model -l 500
16/3: trainer.fit(LightningDMF(model))
16/4: trainer.fit(LightningDMF(model))
16/5: X_tr_shuf.dtype
16/6: X_tr_shuf = X_tr_shuf.float()
16/7: X_val = X_val.float()
16/8: X_tr_shuf.dtype
16/9: trainer.fit(LightningDMF(model))
16/10: Y_tr_shuf.dtype
16/11:
class LightningDMF(pl.LightningModule):
    def __init__(self, model):
        super().__init__()
        self.model = model

    def forward(self, x):
        return self.model(x)

    def training_step(self, batch, batch_idx):
        x, y = batch
        y_hat = self.model(x)
        loss = F.cross_entropy(y_hat, y)
        self.log("train_loss", loss)
        return loss

    def validation_step(self, batch, batch_idx):
        x, y = batch
        y_hat = self.model(x)
        loss = F.cross_entropy(y_hat, y)
        self.log("val_loss", loss)
        return loss

    def configure_optimizers(self):
        return torch.optim.Adam(self.parameters(), lr=1e-3)
16/12: trainer.fit(LightningDMF(model), train_dataloaders=DataLoader(TensorDataset(X_tr_shuf, Y_tr_shuf), batch_size=batch_size, shuffle=True), val_dataloaders=DataLoader(TensorDataset(X_val, Y_val), batch_size=batch_size))
16/13:
from DMF import DMF
trainer.fit(LightningDMF(model), train_dataloaders=DataLoader(TensorDataset(X_tr_shuf, Y_tr_shuf), batch_size=batch_size, shuffle=True), val_dataloaders=DataLoader(TensorDataset(X_val, Y_val), batch_size=batch_size))
16/14:
from DMF import DMF
trainer.fit(LightningDMF(model), train_dataloaders=DataLoader(TensorDataset(X_tr_shuf, Y_tr_shuf), batch_size=batch_size, shuffle=True), val_dataloaders=DataLoader(TensorDataset(X_val, Y_val), batch_size=batch_size))
16/15:
from DMF import DMF
trainer.fit(LightningDMF(model), train_dataloaders=DataLoader(TensorDataset(X_tr_shuf, Y_tr_shuf), batch_size=batch_size, shuffle=True), val_dataloaders=DataLoader(TensorDataset(X_val, Y_val), batch_size=batch_size))
16/16:
from DMF import DMF
trainer.fit(LightningDMF(model), train_dataloaders=DataLoader(TensorDataset(X_tr_shuf, Y_tr_shuf), batch_size=batch_size, shuffle=True), val_dataloaders=DataLoader(TensorDataset(X_val, Y_val), batch_size=batch_size))
16/17:
Y_tr_shuf = Y_tr_shuf.float()
Y_val = Y_val.float()
16/18:
from DMF import DMF
trainer.fit(LightningDMF(model), train_dataloaders=DataLoader(TensorDataset(X_tr_shuf, Y_tr_shuf), batch_size=batch_size, shuffle=True), val_dataloaders=DataLoader(TensorDataset(X_val, Y_val), batch_size=batch_size))
16/19: trainer.fit(LightningDMF(model), train_dataloaders=DataLoader(TensorDataset(X_tr_shuf, Y_tr_shuf), batch_size=batch_size, shuffle=True), val_dataloaders=DataLoader(TensorDataset(X_val, Y_val), batch_size=batch_size))
16/20: trainer.fit(LightningDMF(model), train_dataloaders=DataLoader(TensorDataset(X_tr_shuf, Y_tr_shuf), batch_size=batch_size, shuffle=True), val_dataloaders=DataLoader(TensorDataset(X_val, Y_val), batch_size=batch_size))
17/1: %tensorboard --logdir data/pt_models/lightning_logs
17/2: %load_ext tensorboard
17/3: %tensorboard --logdir data/pt_models/lightning_logs
19/1:
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
from DMF import LightningDMF
model = LightningDMF.load_from_checkpoint('data/pt_logs/checkpoint_prok/epoch=669-step=217080-val_loss=0.03-val_acc=0.98.ckpt', model=DMF_3class(), map_location=device)
19/2:
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
19/3:
import pandas as pd
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
19/4: batch_size = 128
19/5:
import pandas as pd
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
19/6:
target = []
pred = []
for x, y in val_dataloaders:
    target.append(y)
    pred.append(torch.argmax(model(x), dim=1))
19/7: target_1 = torch.concat(target, dim=1)
19/8: target_1 = torch.concat(target)
19/9: pred = torch.concat(pred)
19/10: from torchmetrics.functional import f1_score
19/11: f1_score(target, pred)
19/12: f1_score(target_1, pred)
19/13: f1_score(target_1, pred, average=None)
19/14: f1_score(target_1, pred, average=None, num_classes=3)
19/15: val_dataset.sequence_info[0]
19/16: val_dataset.sequence_info.iloc[0]
19/17: val_dataset.sequence_info[val_dataset.sequence_info['id'] == 'CP035932.1']
19/18: train_dataset.sequence_info[train_dataset.sequence_info['id'] == 'CP035932.1']
19/19: from Bio import SeqIO
19/20: import utils
19/21: records = list(SeqIO.parse('data/sampled_sequence_subsampled/03_PROK450_PROKVIR225_PLSMD225_EUK67_EUKVIR33.fa', 'fasta'))
19/22: r = records[767]
19/23: onehot = utils.seq2onehot(str(r.seq))
19/24: onehot = onehot[:, :4]
19/25: onehot = torch.from_numpy(onehot).float().to(device)
19/26: model(onehot)
19/27: model(onehot[None, None, :, :])
19/28: model(onehot[None, :, :])
19/29: onehot = onehot.to('cpu')
19/30: model(onehot[None, :, :])
19/31: torch.softmax(model(onehot[None, :, :]))
19/32: torch.softmax(model(onehot[None, :, :]), 1)
19/33: onehot = utils.seq2onehot(str(records[773]))
19/34: onehot = utils.seq2onehot(str(records[773]))[:,:4]
19/35: onehot = torch.from_numpy(onehot).float()
19/36: torch.softmax(model(onehot[None, :, :]), 1)
19/37: x = next(val_dataloaders)
19/38: a = torch.softmax(model(x), 1)
19/39: val_dataset.sequence_info.iloc[0]
19/40: val_dataset.sequence_info.iloc[1]
19/41:
import os
os.path.exists('data/single_fasta/NZ_CP041997.1.fasta')
19/42: r = SeqIO.parse('data/single_fasta/NZ_CP041997.1.fasta', 'fasta')
19/43: r = next(r)
19/44: onehot = utils.seq2onehot(r, 5000)
19/45: a = r.seq
19/46: 'N' in a
19/47: onehot = utils.seq2onehot(r, 5000)
19/48: b = np.load('data/single_onehot/NZ_CP041997.1.npy')
19/49: b = np.load('data/single_onehot/NZ_CP041997.1.fasta.npy')
19/50: c = utils.seq2onehot(a)
19/51: (c.sum(axis=1)==0).sum()
19/52: onehot = utils.seq2onehot(r, 5000)
19/53: c = utils.seq2intseq(r, 5000)
19/54: b = utils.int2onehot(c)
19/55: c = utils.seq2intseq(r.seq, 5000)
19/56: onehot = utils.seq2onehot(r.seq, 5000)
20/1: import pandas as pd
20/2: pd.read_csv('data/final_test.csv')
20/3: test = pd.read_csv('data/final_test.csv')
20/4: test['class'].count()
20/5: test['class'].count(axis='columns')
20/6: test['class'].value_counts()
20/7:
import os
for _, row in test.iterrows():
    if not os.path.exists(f'data/single_onehot/{row['id']}.fasta.npy'):
        print(row)
20/8:
import os
for _, row in test.iterrows():
    if not os.path.exists(f'data/single_onehot/{row["id"]}.fasta.npy'):
        print(row)
21/1:
import numpy as np
from matplotlib import pyplot as plt
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
ratio = np.concatenate([np.concatenate([PROK_EUK_RATIO[i, 0] * PROK_PROKVIR_PLASMID_RATIO, PROK_EUK_RATIO[i, 1] * EUK_EUKVIR_RATIO], axis=1) for i in range(PROK_EUK_RATIO.shape[0])], axis=0)
21/2:
import numpy as np
from matplotlib import pyplot as plt
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
ratio = np.concatenate([np.concatenate([PROK_EUK_RATIO[i, 0] * PROK_PROKVIR_PLASMID_RATIO, PROK_EUK_RATIO[i, 1] * EUK_EUKVIR_RATIO], axis=1) for i in range(PROK_EUK_RATIO.shape[0])], axis=0)
21/3:
import pandas as pd
ratio_df = pd.DataFrame(ratio, columns=['Prok', 'ProkVirus', 'Plasmid', 'Euk', 'EukVirus'])
ratio_df.plot(kind='bar', stacked=True, ax=ax, color=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd'])
21/4:
fig, ax = plt.subplots()

ratio_df = pd.DataFrame(ratio, columns=['Prok', 'ProkVirus', 'Plasmid', 'Euk', 'EukVirus'])
ratio_df.plot(kind='bar', stacked=True, ax=ax, color=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd'])
21/5: ratio_df.plot(kind='bar', stacked=True, ax=ax)
21/6: ratio_df.plot(kind='bar', stacked=True)
21/7:
ratio_df.plot(kind='bar',
    stacked=True,
    # ax=ax,
    width=0.8,
    xaxis='off',
    color=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    )
21/8:
ratio_df.plot(kind='bar',
    stacked=True,
    # ax=ax,
    width=0.8,
    xticks=None,
    color=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    )
21/9:
ratio_df.plot(kind='bar',
    stacked=True,
    # ax=ax,
    width=0.8,
    xticks=None,
    color=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    )
plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
21/10:
ratio_df.plot(kind='bar',
    stacked=True,
    # ax=ax,
    width=0.8,
    xticks=None,
    color=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    )
plt.legend(loc='upper left')
21/11:
ratio_df.plot(kind='bar',
    stacked=True,
    # ax=ax,
    width=0.8,
    xticks=None,
    color=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    )
plt.legend(loc='upper left', ncol=5)
21/12:
ratio_df.plot(kind='bar',
    stacked=True,
    # ax=ax,
    width=0.8,
    xticks=None,
    color=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    )
plt.legend(loc='upper left', ncol=5, bbox_to_anchor=(0, 1.02, 1, 0.2), mode='expand', borderaxespad=0, frameon=False)
21/13:
ratio_df.plot(kind='bar',
    stacked=True,
    # ax=ax,
    width=0.8,
    xticks=None,
    color=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    )
plt.legend(loc='upper center', ncol=5, bbox_to_anchor=(0, 1.02, 1, 0.2), mode='expand', borderaxespad=0, frameon=False)
21/14: plt.legend(loc='upper center', ncol=5, bbox_to_anchor=(0, 1, 1, 0.2), mode='expand', borderaxespad=0, frameon=False)
21/15: plt.legend(loc='upper center', ncol=5, bbox_to_anchor=(0, 0.5, 1, 0.2), mode='expand', borderaxespad=0, frameon=False)
21/16:
ratio_df = pd.DataFrame(ratio, columns=['Prok', 'ProkVirus', 'Plasmid', 'Euk', 'EukVirus'])
ratio_df.plot(kind='bar',
    stacked=True,
    # ax=ax,
    width=0.8,
    xticks=None,
    color=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    )
plt.legend(loc='upper center', ncol=5, bbox_to_anchor=(0, 0.5, 1, 0.2), mode='expand', borderaxespad=0, frameon=False)
21/17:
ratio_df.plot(kind='bar',
    stacked=True,
    # ax=ax,
    width=0.8,
    xticks=None,
    color=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    )
plt.legend(loc='upper center', ncol=5, bbox_to_anchor=(0, 1, 0.8, 0.2), mode='expand', borderaxespad=0, frameon=False)
21/18:
ratio_df = pd.DataFrame(ratio, columns=['Prok', 'ProkVirus', 'Plasmid', 'Euk', 'EukVirus'])
ratio_df.plot(kind='bar',
    stacked=True,
    # ax=ax,
    width=0.8,
    xticks=None,
    color=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    )
plt.legend(loc='upper center', ncol=5, bbox_to_anchor=(0, 1, 1, 0.2), borderaxespad=0, frameon=False)
21/19:
ratio_df.plot(kind='bar',
    stacked=True,
    # ax=ax,
    width=0.8,
    xticks=None,
    color=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    )
plt.legend(loc='upper center', ncol=5, bbox_to_anchor=(0, 1, 1, 0.2), mode='expand', borderaxespad=0, frameon=False)
21/20:
ratio_df.plot(kind='bar',
    stacked=True,
    # ax=ax,
    width=0.8,
    xticks=None,
    color=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    )

plt.legend(loc='upper center', ncol=5, bbox_to_anchor=(0, 1, 1, 0.2), mode='expand', borderaxespad=0, frameon=True)
21/21:
ratio_df.plot(kind='bar',
    stacked=True,
    # ax=ax,
    width=0.8,
    xticks=None,
    color=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    )

plt.legend(loc='upper center', ncol=5, bbox_to_anchor=(0, 1, 1, 0.2), mode='expand', borderaxespad=0, frameon=True, labelspacing=0.1,)
21/22:
ratio_df.plot(kind='bar',
    stacked=True,
    # ax=ax,
    width=0.8,
    xticks=None,
    color=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    )
plt.legend(loc='upper center', ncol=5, bbox_to_anchor=(0, 1, 1, 0.2), mode='expand', borderaxespad=0, frameon=True, labelspacing=0.1, columnspacing=0.5, handletextpad=0.2, handlelength=1.5)
21/23: acc_df = pd.DataFrame([tiara_acc, whok_acc, dmf_acc], columns=['Tiara', 'WHO-K', 'DMF'])
21/24:
tiara_acc = [0.989010989010989, 0.989, 0.989, 0.993, 0.979, 0.978021978021978, 0.984, 0.988, 0.9669669669669669, 0.964964964964965, 0.966, 0.979, 0.953, 0.966, 0.962, 0.952, 0.92992992992993, 0.945054945054945, 0.945, 0.941]
tiara_f1 = [0.9299363057324841, 0.9271523178807948, 0.9219858156028369, 0.9448818897637796, 0.9561586638830899, 0.9521739130434782, 0.9634703196347033, 0.9690721649484536, 0.9588014981273408, 0.9542483660130718, 0.9527777777777778, 0.9674418604651163, 0.9580731489741303, 0.96875, 0.9624505928853755, 0.945945945945946, 0.9510489510489509, 0.9602888086642599, 0.9575289575289575, 0.9483814523184603,]

whok_acc = [0.9585062240663901, 0.9611344537815126, 0.9512448132780082, 0.9473140495867769, 0.9332591768631813, 0.9394618834080718, 0.9144444444444444, 0.9058956916099773, 0.9076354679802956, 0.882640586797066, 0.8824969400244798, 0.8544776119402985, 0.8688741721854305, 0.8481182795698925, 0.8297872340425532, 0.789261744966443, 0.8418740849194729, 0.8179148311306902, 0.7800586510263929, 0.7246376811594203,]
whok_f1 = [0.9777777777777779, 0.9793411501954216, 0.9737576772752652, 0.9716509171762089, 0.9582753824756607, 0.9626038781163434, 0.9466389466389465, 0.9439567859554355, 0.9298409728718428, 0.9114391143911438, 0.9136690647482014, 0.8963684676705048, 0.8626907073509016, 0.844566712517194, 0.8302387267904509, 0.8110709987966306, 0.686046511627907, 0.6787564766839378, 0.6411483253588516, 0.6169354838709677,]

dmf_acc = [0.998001998001998, 0.998, 1.0, 0.998, 0.996, 0.999000999000999, 0.998, 0.998, 0.995995995995996, 0.995995995995996, 0.995, 0.995, 0.994, 0.991, 0.996, 0.991, 0.992992992992993, 0.993006993006993, 0.991, 0.987,]
dmf_f1 = [0.9878048780487805, 0.9876543209876543, 1.0, 0.9850746268656716, 0.992, 0.997920997920998, 0.9955357142857144, 0.995, 0.9951923076923077, 0.9949874686716792, 0.9933065595716198, 0.992503748125937, 0.9948453608247423, 0.9919282511210762, 0.9961977186311786, 0.9903743315508022, 0.9953364423717521, 0.9951287404314543, 0.9933184855233853, 0.9891213389121339,]
acc_df = pd.DataFrame([tiara_acc, whok_acc, dmf_acc], columns=['Tiara', 'WHO-K', 'DMF'])
21/25:
tiara_acc = np.array(tiara_acc)
tiara_f1 = np.array(tiara_f1)
whok_acc = np.array(whok_acc)
whok_f1 = np.array(whok_f1)
dmf_acc = np.array(dmf_acc)
dmf_f1 = np.array(dmf_f1)
21/26: acc_df = pd.DataFrame([tiara_acc, whok_acc, dmf_acc], columns=['Tiara', 'WHO-K', 'DMF'])
21/27: acc_df = pd.DataFrame(columns=['Tiara', 'WHO-K', 'DMF'])
21/28:
acc_df['Tiara'] = tiara_acc
acc_df['WHO-K'] = whok_acc
acc_df['DMF'] = dmf_acc
21/29: acc_df.plot(kind='bar')
21/30:
acc_df.plot(kind='bar',
    width=0.8,
)
21/31:
acc_df.plot(kind='bar',
    width=0.8,
)
plt.legend(loc='upper left')
21/32:
ratio_df = pd.DataFrame(ratio, columns=['Prok', 'ProkVirus', 'Plasmid', 'Euk', 'EukVirus'])
ratio_df.plot(kind='bar',
    stacked=True,
    # ax=ax,
    width=0.8,
    xticks=None,
    color=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    )
21/33:
xticks = [f'DS_{i}' for i in range(1, 21)]
ratio_df.plot(kind='bar',
    stacked=True,
    # ax=ax,
    width=0.8,
    xticks=xla,
    color=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    )
21/34:
xticks = [f'DS_{i}' for i in range(1, 21)]
ratio_df.plot(kind='bar',
    stacked=True,
    # ax=ax,
    width=0.8,
    xticks=xticks,
    color=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    )
21/35:
ratio_df = pd.DataFrame(ratio, columns=['Prok', 'ProkVirus', 'Plasmid', 'Euk', 'EukVirus'], index=xticks)
ratio_df.plot(kind='bar',
    stacked=True,
    # ax=ax,
    width=0.8,
    # xticks=xla,
    color=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    )
21/36:
ratio_df.plot(kind='bar',
    stacked=True,
    # ax=ax,
    width=0.8,
    # xticks=xla,
    color=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    )
plt.xticks(rotation=45)
21/37:
ratio_df.plot(kind='bar',
    stacked=True,
    # ax=ax,
    width=0.8,
    # xticks=xla,
    color=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    )
plt.xticks(rotation=45)
plt.margins(0.01)
21/38:
acc_df = pd.DataFrame(columns=['Tiara', 'Whok', 'DMF'])
acc_df['Tiara'] = tiara_acc
acc_df['Whok'] = whok_acc
acc_df['DMF'] = dmf_acc

acc_df.plot(kind='bar',
    width=0.8,
)
plt.legend(loc='upper left')
21/39:
acc_df = pd.DataFrame(columns=['Tiara', 'Whok', 'DMF'])
acc_df['Tiara'] = tiara_acc
acc_df['Whok'] = whok_acc
acc_df['DMF'] = dmf_acc

acc_df.plot(kind='bar',
    width=0.8,
)
plt.legend(loc='lower left')
plt.grid(axis='y')
21/40:
acc_df = pd.DataFrame(columns=['Tiara', 'Whok', 'DMF'])
acc_df['Tiara'] = tiara_acc
acc_df['Whok'] = whok_acc
acc_df['DMF'] = dmf_acc

ax = acc_df.plot(kind='bar',
    width=0.8,
)
ax.legend(loc='lower left')
ax.grid(axis='y')
21/41:
acc_df = pd.DataFrame(columns=['Tiara', 'Whok', 'DMF'])
acc_df['Tiara'] = tiara_acc
acc_df['Whok'] = whok_acc
acc_df['DMF'] = dmf_acc

ax = acc_df.plot(kind='bar',
    width=0.8,
)
ax.legend(loc='lower left')
ax.grid(axis='y')
ax.set_axisbelow(True)
21/42:
acc_df = pd.DataFrame(columns=['Tiara', 'Whok', 'DMF'])
acc_df['Tiara'] = tiara_acc
acc_df['Whok'] = whok_acc
acc_df['DMF'] = dmf_acc

ax = acc_df.plot(kind='bar',
    width=0.8,
)
ax.legend(loc='lower left')
ax.grid(axis='y')
ax.set_axisbelow(True)
ax.margins(0.01)
21/43:
acc_df = pd.DataFrame(columns=['Tiara', 'Whok', 'DMF'])
acc_df['Tiara'] = tiara_acc
acc_df['Whok'] = whok_acc
acc_df['DMF'] = dmf_acc

ax = acc_df.plot(kind='bar',
    width=0.8,
)
ax.legend(loc='lower left')
ax.grid(axis='y')
ax.set_axisbelow(True)
ax.margins(0.0)
21/44:
import seaborn as sns
sns.catplot(data=acc_df, kind='bar',
    height=4,
    aspect=2
    )
21/45:
ratio_df = pd.DataFrame(ratio, columns=['Prok', 'ProkVirus', 'Plasmid', 'Euk', 'EukVirus'], index=xticks)
ratio_df.plot(kind='bar',
    stacked=True,
    # ax=ax,
    width=0.8,
    # xticks=xla,
    color=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    )
plt.xticks(rotation=45)
ax.grid(axis='y')
ax.set_axisbelow(True)
plt.margins(0.01)
21/46:
RESULT_DIR = 'data/result_plasflow'
results_fn = os.listdir(RESULT_DIR)
results_fn = [f'{f}/whokaryote_predictions_T.tsv' for f in results_fn if not f.startswith('log')]
results_fn.sort()
result = pd.read_table(os.path.join(RESULT_DIR, next(results_fn)), index_col=0)
21/47:
RESULT_DIR = 'data/result_plasflow'
results_fn = os.listdir(RESULT_DIR)
results_fn = [f'{f}/whokaryote_predictions_T.tsv' for f in results_fn if not f.startswith('log')]
results_fn.sort()
result = pd.read_table(os.path.join(RESULT_DIR, results_fn[0]), index_col=0)
21/48:
RESULT_DIR = 'data/result_plasflow'
results_fn = os.listdir(RESULT_DIR)
results_fn.sort()
result = pd.read_table(os.path.join(RESULT_DIR, results_fn[0]), index_col=0)
21/49:
RESULT_DIR = 'data/result_plasflow'
results_fn = os.listdir(RESULT_DIR)
results_fn.sort()
result = pd.read_table(os.path.join(RESULT_DIR, results_fn[0]), index_col=2)
21/50: result.iloc[28,:]['label']
21/51: result.iloc[28,:]['label'].startswith('plasmid')
21/52: result = pd.read_table('data/result_ppr/00_PROK643_PROKVIR129_PLSMD129_EUK83_EUKVIR17.csv', sep=',')
21/53:
index = list(result['Header'])
index = [i[1:].split()[0] for i in index]
result.index = index
21/54: result = pd.read_table('data/result_ppr/00_PROK643_PROKVIR129_PLSMD129_EUK83_EUKVIR17.csv', sep=',')
21/55:
index = list(result['Header'])
index = [i[1:].split()[0] for i in index]
result.index = index
21/56: result = pd.read_table('data/result_ppr/00_PROK643_PROKVIR129_PLSMD129_EUK83_EUKVIR17.csv', sep=',')
21/57:
index = list(result['Header'])
index = [isplit()[0] for i in index]
result.index = index
21/58:
index = list(result['Header'])
index = [i.split()[0] for i in index]
result.index = index
21/59: len('.fa_gt1bp_dvfpred.txt')
21/60: result = pd.read_table('data/result_virsorter2/00_PROK643_PROKVIR129_PLSMD129_EUK83_EUKVIR17/final-viral-score.tsv')
21/61: result = pd.read_table('data/result_virsorter2/00_PROK643_PROKVIR129_PLSMD129_EUK83_EUKVIR17/final-viral-score.tsv', index_col=0)
21/62:
index = list(result.index)
    index = [i.split('|')[0] for i in index]
    result.index = index
21/63:
index = list(result.index)
index = [i.split('|')[0] for i in index]
result.index = index
21/64: len('VIBRANT_phages_')
21/65: result = pd.read_table('data/result_vibrant/VIBRANT_phages_00_PROK643_PROKVIR129_PLSMD129_EUK83_EUKVIR17/00_PROK643_PROKVIR129_PLSMD129_EUK83_EUKVIR17.phages_combined.txt', index_col=0, header=None)
21/66:
index = list(result.index)
index = [i.split()[0] for i in index]
result.index = index
21/67: 'VIBRANT_phages_00_PROK643_PROKVIR129_PLSMD129_EUK83_EUKVIR17'[16:]
21/68: %run scripts/plot_auroc.py
21/69: %run scripts/plot_auroc.py
21/70: %run scripts/plot_auroc.py
21/71: %run scripts/plot_auroc.py
21/72: %run scripts/plot_auroc.py
21/73: %run scripts/plot_auroc.py
21/74: %run scripts/plot_auroc.py
21/75: %run scripts/plot_auroc.py
21/76: %run scripts/plot_auroc.py
21/77: %run scripts/plot_auroc.py
21/78: %run scripts/plot_auroc.py
21/79: %run scripts/plot_auroc.py
21/80: %run scripts/plot_auroc.py
21/81: %run scripts/plot_auroc.py
21/82: %run scripts/plot_auroc.py
21/83: %run scripts/plot_auroc.py
21/84: %run scripts/plot_auroc.py
21/85: %run scripts/plot_auroc.py
21/86: %run scripts/plot_auroc.py
21/87: %run scripts/plot_auroc.py
21/88:

%run scripts/plot_auroc.py
21/89:

%run scripts/plot_auroc.py
21/90: %run scripts/plot_auroc.py
21/91:

%run scripts/plot_auroc.py
21/92: %run scripts/plot_auroc.py
21/93: %run scripts/plot_auroc.py
21/94: %run scripts/plot_auroc.py
21/95: %run scripts/plot_auroc.py
21/96:
from Bio import SeqIO
from tqdm import tqdm
a = [len(record) for record in tqdm(SeqIO.parse('data/SPOT/SPOT_5YrDuraVir_min_5000_newbler_toAmos_minimus2_id0.99_renamed.fa', 'fasta'))]
21/97: a = np.array(a)
21/98: a.mean()
21/99: a.std()
21/100: np.median(a)
21/101: plasmid_meta = pd.read_table('data/plsdb/plsdb.tsv')
21/102: plasmid_meta.columns
21/103:
b = [len(record) for record in tqdm(SeqIO.parse('/home/shengwei/Viromes/04_read_mapping_per_sample/assemblies_per_sample/2SPOviralmetaG_derep.fasta', 'fasta'))]
b = np.array(b)
21/104: b.mean()
21/105: np.median(b)
21/106: spot_output = pd.read_table('data/SPOT/SPOT_5YrDuraVir_min_5000_newbler_toAmos_minimus2_id0.99_renamed_pred.csv', sep=',')
21/107: spot_class = spot_output.argmax(axis=1)
21/108: spot_class = spot_output.iloc[:, 1:].to_numpy()
21/109: spot_class = spot_class.argmax(axis=1)
21/110: proportion = [sum(spot_class==i)/len(spot_class) for i in range(5)]
21/111: proportion
21/112: sum(proportion)
21/113: pd.read_csv('data/final_test.csv')
21/114:
import pickle
data = pickle.load('results/contig_len_50000.pkl')
21/115:
import pickle
data = pickle.load(open('results/contig_len_50000.pkl', 'rb'))
21/116:
all_y = data['all_y']
all_y_hat = data['all_y_hat']
21/117: all_y_hat = all_y_hat.reshape(-1, 5)
21/118: print(auroc(torch.tensor(all_y_hat), torch.tensor(all_y), num_classes=5, average='None'))
21/119: %run scripts/plot_auroc.py
21/120: pickle.dump({'all_y': all_y, 'all_y_hat': all_y_hat}, open('results/contig_len_50000.pkl', 'wb'))
21/121: pickle.dump({'all_y': all_y, 'all_y_hat': all_y_hat}, open('results/contig_len_50000.pkl', 'wb'))
21/122: %run scripts/plot_auroc.py
21/123: fpr, tpr, _ = roc(torch.tensor(all_y_hat), torch.tensor(all_y), num_classes=5)
21/124: data = pickle.load(open('results/contig_len_50000.pkl', 'rb'))
21/125: all_y_hat = all_y_hat.reshape(-1, 5)
21/126: pickle.dump({'all_y': all_y, 'all_y_hat': all_y_hat}, open('results/contig_len_50000.pkl', 'wb'))
21/127:
data = pickle.load(open('results/contig_len_50000.pkl', 'rb'))
all_y = data['all_y']
all_y_hat = data['all_y_hat']
21/128: %run scripts/plot_auroc.py
21/129: del a all_y all_y_hat b fpr tpr
21/130: del a,all_y,all_y_hat, b, fpr, tpr
21/131: test = pd.read_csv('data/final_test.csv')
21/132: test['class'].count()
21/133: np.array(test['class']).unique()
21/134: np.unique(np.array(test['class']))
21/135: np.unique(np.array(test['class']), return_counts=True)
21/136: _, a = np.unique(np.array(test['class']), return_counts=True)
21/137: a
21/138: a.T
21/139: np.unique(all_y, return_counts=True)
21/140:
data = pickle.load(open('results/contig_len_50000.pkl', 'rb'))
all_y = data['all_y']
all_y_hat = data['all_y_hat']
21/141: np.unique(all_y, return_counts=True)
21/142: _, spot_count = np.unique(spot_output, return_counts=True)
21/143: spot_output.max()
21/144: _, spot_count = np.unique(spot_class, return_counts=True)
21/145: spot_count / spot_count.sum()
21/146: results = [pd.read_csv('results/results_mu0_delta0.csv'), pd.read_csv('results/results_mu0.02_delta0.02.csv'), pd.read_csv('results/results_mu0.05_delta0.05.csv')]
21/147: results = [pd.read_csv('results/results_mu0_delta0.csv'), pd.read_csv('results/results_mu0.002_delta0.002.csv'), pd.read_csv('results/results_mu0.005_delta0.005.csv')]
21/148: results = [pd.read_csv('results/results_mu0_delta0.csv'), pd.read_csv('results/results_mu0_delta0.05.csv'), pd.read_csv('results/results_mu0.05_delta0.csv'), pd.read_csv('results/results_mu0.005_delta0.005.csv')]
21/149:
out = np.concatenate([result.iloc[:, i] for result in results])
df = pd.DataFrame(out, columns=['mu=0, delta=0', 'mu=0, delta=0.05', 'mu=0.05, delta=0', 'mu=0.05, delta=0.05'], index=xticks)
21/150:
out = np.concatenate([result.iloc[:, 0] for result in results])
df = pd.DataFrame(out, columns=['mu=0, delta=0', 'mu=0, delta=0.05', 'mu=0.05, delta=0', 'mu=0.05, delta=0.05'], index=xticks)
21/151:
out = np.concatenate([result.iloc[:, i] for result in results], axis=1)
df = pd.DataFrame(out, columns=['mu=0, delta=0', 'mu=0, delta=0.05', 'mu=0.05, delta=0', 'mu=0.05, delta=0.05'], index=xticks)
21/152:
out = np.concatenate([result.iloc[:, 0] for result in results], axis=1)
df = pd.DataFrame(out, columns=['mu=0, delta=0', 'mu=0, delta=0.05', 'mu=0.05, delta=0', 'mu=0.05, delta=0.05'], index=xticks)
21/153: b = [result.iloc[:, 0] for result in results]
21/154: b[0].shape
21/155: out = np.concatenate([result.iloc[:, 0][None, :] for result in results], axis=1)
21/156: out = np.concatenate([result.iloc[:, 0].to_numpy()[:, None] for result in results], axis=1)
21/157: df = pd.DataFrame(out, columns=['mu=0, delta=0', 'mu=0, delta=0.05', 'mu=0.05, delta=0', 'mu=0.05, delta=0.05'], index=xticks)
21/158: %run scripts/plot_mutation.py
21/159: %run scripts/plot_mutation.py
21/160: %run scripts/plot_mutation.py
21/161: %run scripts/plot_mutation.py
21/162: %run scripts/plot_mutation.py
21/163: %run scripts/plot_mutation.py
   1:
import pandas as pd
test = pd.read_csv('data/final_test.csv')
train = pd.read_csv('data/train.csv')
val = pd.read_csv('data/test.csv')
   2: test.index.intersection(train.index)
   3: len(test.index.intersection(train.index))
   4: len(test.index.intersection(val.index))
   5: len(train.index.intersection(val.index))
   6: len(set(test.index).intersection(train.index))
   7: len(set(val.index).intersection(train.index))
   8: set(val.index)a =
   9: a = set(val.index)
  10:
test = pd.read_csv('data/final_test.csv', index_col=0)
train = pd.read_csv('data/train.csv', index_col=0)
val = pd.read_csv('data/test.csv', index_col=0)
  11: a = set(val.index)
  12: test.index.intersection(train.index)
  13: train.index.intersection(val.index)
  14: test.index.intersection(val.index)
  15: %hist -o -g -f ipython_history.md
