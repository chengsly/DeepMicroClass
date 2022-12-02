import itertools
import torch
import torch.nn as nn
import torch.nn.functional as F
import math
import pytorch_lightning as pl
from torchmetrics.functional import accuracy, f1_score, auroc
from torch import Tensor

# Channel configuration
# Res layer number

class LightningDMF(pl.LightningModule):
    def __init__(self, model, lr=1e-3, weight_decay=1e-5, weight=None, num_classes=5, **kwargs):
        super().__init__()
        self.model = model
        # self.weight = weight
        self.num_classes = num_classes

        self.lr = lr
        self.weight_decay = weight_decay

        self.save_hyperparameters(ignore=['model'])

    def forward(self, x):
        return self.model(x)

    def training_step(self, batch, batch_idx):
        x, y = batch
        y_hat = self.model(x)
        loss = F.cross_entropy(y_hat, y)
        self.log("train_loss", loss)
        self.log('train_acc_epoch', accuracy(y_hat, y.int(), average='macro', num_classes=self.num_classes), on_step=False, on_epoch=True, prog_bar=True)
        self.log('train_f1_epoch', f1_score(y_hat, y.int(), average='macro', num_classes=self.num_classes), on_step=False, on_epoch=True, prog_bar=True)
        return loss

    def validation_step(self, batch, batch_idx):
        x, y = batch
        y_hat = self.model(x)
        val_loss = F.cross_entropy(y_hat, y)
        self.log("val_loss", val_loss, prog_bar=True)
        self.log('val_acc', accuracy(y_hat, y.int(), average='macro', num_classes=self.num_classes), prog_bar=True)
        self.log('val_f1', f1_score(y_hat, y.int(), average='macro', num_classes=self.num_classes), prog_bar=True)
        return val_loss
    
    def test_step(self, batch, batch_idx):
        loss, acc, auc = self._shared_eval_step(batch, batch_idx)
        metrics = {"test_acc": acc, "test_loss": loss}
        self.log_dict(metrics)
        return metrics

    def _shared_eval_step(self, batch, batch_idx):
        x, y = batch
        y_hat = self.model(x)
        loss = F.cross_entropy(y_hat, y)
        acc = accuracy(y_hat, y)
        auc = auroc(y_hat, y, num_classes=self.num_classes, average=None)
        return loss, acc, auc

    def configure_optimizers(self):
        return torch.optim.Adam(self.parameters(), lr=self.lr, weight_decay=self.weight_decay)


class DeepMicroClass(nn.Module):
    def __init__(self) -> None:
        super().__init__()

        self.codon_transformer = CodonTransformer()

        self.fc = nn.Sequential(
            nn.Linear(1024, 256),
            # nn.Linear(512, 256),
            nn.PReLU(),
            nn.Dropout(p=0.2, inplace=True),
            nn.Linear(256, 5)
        )

        self.base_channel = nn.Sequential(
            nn.Conv2d(1, 64, (6, 4)),
            nn.PReLU(),
            nn.Flatten(start_dim=2),
            nn.AvgPool1d(3),
            nn.BatchNorm1d(64),
            nn.Conv1d(64, 128, 3),
            nn.PReLU(),
            nn.AvgPool1d(3),
            nn.BatchNorm1d(128),
            nn.Conv1d(128, 256, 3),
            nn.PReLU(),
            nn.BatchNorm1d(256),
            nn.AdaptiveAvgPool1d(1),
            nn.Flatten(),
        )

        self.codon_channel = nn.Sequential(
            nn.Conv2d(1, 64, (2, 64)),
            nn.PReLU(),
            nn.Flatten(start_dim=2),
            # nn.MaxPool1d(3),
            nn.AvgPool1d(3),
            nn.BatchNorm1d(64),
            nn.Conv1d(64, 128, 3),
            nn.PReLU(),
            # nn.MaxPool1d(3),
            nn.AvgPool1d(3),
            nn.BatchNorm1d(128),
            nn.Conv1d(128, 256, 3),
            nn.PReLU(),
            nn.BatchNorm1d(256),
            # nn.Conv1d(256, 5, 1),
            nn.AdaptiveAvgPool1d(1),
            nn.Flatten(),
        )

        self.dropout = nn.Dropout(p=0.2, inplace=True)

    def forward(self, x):
        interm_forward = self.base_channel(x)
        codon = self.codon_transformer(x)
        codon = self.codon_channel(codon)
        rev = torch.flip(x, dims=[-1, -2])
        interm_backward = self.base_channel(rev)
        codon_backward = self.codon_transformer(rev)
        codon_backward = self.codon_channel(codon_backward)
        z = torch.cat((interm_forward, codon, interm_backward, codon_backward, ), dim=1)
        # z = torch.cat((codon, codon_backward, ), dim=1)
        # z = torch.cat((interm_forward, interm_backward, ), dim=1)
        z = self.dropout(z)
        z = self.fc(z)
        # z = (codon + codon_backward) / 2
        # z = codon
        return z

class DeepMicroClass_base(nn.Module):
    def __init__(self, lr=1e-3, k=6) -> None:
        super().__init__()

        self.lr = lr
        self.k = k
        
        self.fc = nn.Sequential(
            nn.Linear(512, 256),
            nn.PReLU(),
            nn.Dropout(p=0.2, inplace=True),
            nn.Linear(256, 5)
        )

        self.base_channel = nn.Sequential(
            nn.Conv2d(1, 64, (6, 4)),
            nn.PReLU(),
            nn.Flatten(start_dim=2),
            nn.AvgPool1d(3),
            nn.BatchNorm1d(64),
            nn.Conv1d(64, 128, 3),
            nn.PReLU(),
            nn.AvgPool1d(3),
            nn.BatchNorm1d(128),
            nn.Conv1d(128, 256, 3),
            nn.PReLU(),
            nn.BatchNorm1d(256),
            nn.AdaptiveAvgPool1d(1),
            nn.Flatten(),
        )

    def forward(self, x):
        interm_forward = self.base_channel(x)
        rev = torch.flip(x, dims=[-1, -2])
        interm_backward = self.base_channel(rev)
        z = torch.cat((interm_forward, interm_backward, ), dim=1)
        z = self.fc(z)
        return z

class DeepMicroClass_kmer(nn.Module):
    def __init__(self, lr=1e-3, k=3) -> None:
        super().__init__()
        
        self.lr = lr
        self.k = k

        self.fc = nn.Sequential(
            nn.Linear(512, 256),
            nn.PReLU(),
            nn.Dropout(p=0.2, inplace=True),
            nn.Linear(256, 5)
        )

        self.kmer_channel = nn.Sequential(
            nn.Conv2d(1, 64, (2, 64)),
            nn.PReLU(),
            nn.Flatten(start_dim=2),
            nn.AvgPool1d(3),
            nn.BatchNorm1d(64),
            nn.Conv1d(64, 128, 3),
            nn.PReLU(),
            nn.AvgPool1d(3),
            nn.BatchNorm1d(128),
            nn.Conv1d(128, 256, 3),
            nn.PReLU(),
            nn.BatchNorm1d(256),
            nn.AdaptiveAvgPool1d(1),
            nn.Flatten(),
        )

    def forward(self, x):
        kmer = self.kmer_channel(x)
        rev = torch.flip(x, dims=[-1, -2])
        kmer_backward = self.kmer_channel(rev)
        z = torch.cat((kmer, kmer_backward, ), dim=1)
        z = self.fc(z)
        return z


class CodonTransformer(nn.Module):
    def __init__(self) -> None:
        super().__init__()
        self.codon_channel_num = 64
        self.codon_transformer = torch.zeros(64, 1, 3, 4)
        indicies = itertools.product(range(4), repeat=3)
        for i in range(self.codon_transformer.shape[0]):
            index = next(indicies)
            for j in range(3):
                self.codon_transformer[i, 0, j, index[j]] = 1
        self.codon_transformer = nn.Parameter(self.codon_transformer, requires_grad=False)
        self.padding_layers = [nn.ZeroPad2d((0, 0, 0, 2)), nn.ZeroPad2d((0, 0, 0, 1))]


    def forward(self, x):
        mod_len = int(x.shape[2] % 3)
        if mod_len != 2:
            x = self.padding_layers[mod_len](x)
        x = F.conv2d(x, self.codon_transformer) - 2
        x = F.relu(x)
        x = x.flatten(start_dim=2)

        x = x.view(-1, self.codon_channel_num, int(x.shape[2]//3), 3)
        x = x.transpose(2, 3)
        x = x.reshape(-1, 1, self.codon_channel_num, x.shape[-1]*3).transpose(2, 3)
        return x

class KMerTransformer(nn.Module):
    def __init__(self, k=3, rearrange=False) -> None:
        super().__init__()
        self.k = k
        self.rearrange = rearrange
        self.kmer_channel_num = 4 ** self.k
        self.kmer_transformer = torch.zeros(self.kmer_channel_num, 1, self.k, 4)
        indicies = itertools.product(range(4), repeat=self.k)
        for i in range(self.kmer_transformer.shape[0]):
            index = next(indicies)
            for j in range(self.k):
                self.kmer_transformer[i, 0, j, index[j]] = 1
        self.kmer_transformer = nn.Parameter(self.kmer_transformer, requires_grad=False)
        self.padding_layers = [nn.ZeroPad2d((0, 0, 0, i)) for i in range(self.k-1, 0, -1)]


    def forward(self, x):
        mod_len = int(x.shape[2] % self.k)
        if mod_len != self.k-1:
            x = self.padding_layers[mod_len](x)
        x = F.conv2d(x, self.kmer_transformer) - 2
        x = F.relu(x)
        x = x.flatten(start_dim=2)
        if self.rearrange:
            x = x.view(-1, self.kmer_channel_num, int(x.shape[2]//self.k), self.k)
            x = x.transpose(2, 3)
            x = x.reshape(-1, 1, self.kmer_channel_num, x.shape[-1]*self.k).transpose(2, 3)
        return x

class DMFTransformer(nn.Module):
    def __init__(self) -> None:
        super().__init__()

        self.k = 6
        self.kmer_transformer = KMerTransformer(k=self.k, rearrange=True)
        self.embed_d = 128
        self.embedding = nn.Embedding(4 ** self.k, self.embed_d)
        self.pos_encoder = PositionalEncoding(self.embed_d, 0.1)

        self.encoder_layer = nn.TransformerEncoderLayer(d_model=self.embed_d, nhead=4, dim_feedforward=512, dropout=0.1)
        self.transformer_encoder = nn.TransformerEncoder(self.encoder_layer, num_layers=2)
        self.decoder = nn.Linear(self.embed_d, 5)

        self.init_weights()

    def init_weights(self) -> None:
        initrange = 0.1
        self.embedding.weight.data.uniform_(-initrange, initrange)
        self.decoder.bias.data.zero_()
        self.decoder.weight.data.uniform_(-initrange, initrange)
    
    def forward(self, x):
        x = self.kmer_transformer(x,)
        x = torch.argmax(x, dim=1).long()
        x = self.embedding(x) * math.sqrt(self.embed_d)
        x = self.pos_encoder(x)
        z = self.transformer_encoder(x)
        z = z.mean(dim=1)
        z = self.decoder(z)
        return z

def generate_square_subsequent_mask(sz: int) -> Tensor:
    """Generates an upper-triangular matrix of -inf, with zeros on diag."""
    return torch.triu(torch.ones(sz, sz) * float('-inf'), diagonal=1)

class PositionalEncoding(nn.Module):

    def __init__(self, d_model: int, dropout: float = 0.1, max_len: int = 5000):
        super().__init__()
        self.dropout = nn.Dropout(p=dropout)

        position = torch.arange(max_len).unsqueeze(1)
        div_term = torch.exp(torch.arange(0, d_model, 2) * (-math.log(10000.0) / d_model))
        pe = torch.zeros(max_len, 1, d_model)
        pe[:, 0, 0::2] = torch.sin(position * div_term)
        pe[:, 0, 1::2] = torch.cos(position * div_term)
        self.register_buffer('pe', pe)

    def forward(self, x: Tensor) -> Tensor:
        """
        Args:
            x: Tensor, shape [seq_len, batch_size, embedding_dim]
        """
        x = x + self.pe[:x.size(0)]
        return self.dropout(x)

class DMCLSTM(nn.Module):
    def __init__(self, k=3) -> None:
        super().__init__()
        self.k = k
        self.kmer_transformer = KMerTransformer(k=self.k, rearrange=False)
        self.embed_d = 128
        self.embedding = nn.Embedding(4 ** self.k, self.embed_d)
        self.lstm_hidden = 128
        self.lstm = nn.LSTM(self.embed_d, self.lstm_hidden, num_layers=2, batch_first=True, bidirectional=True)
        self.fc = nn.Linear(256, 5)
        self.dropout = nn.Dropout(p=0.5, inplace=True)

    def forward(self, x):
        x = self.kmer_transformer(x)
        x = torch.argmax(x, dim=1).long()
        x = self.embedding(x)# * math.sqrt(self.embed_d)
        _, (x, _) = self.lstm(x)
        # x = x.view(x.size(0), -1)
        x = torch.cat((x[-2], x[-1]), dim=1)
        x = self.dropout(x)
        x = self.fc(x)
        return x


class DMF_tfidf(nn.Module):
    def __init__(self) -> None:
        super().__init__()

        self.codon_transformer = CodonTransformer()

        self.fc = nn.Sequential(
            nn.LazyLinear(512),
            nn.ReLU(inplace=True),
            nn.Dropout(p=0.5),
            nn.Linear(512, 256),
            nn.ReLU(inplace=True),
            nn.Dropout(p=0.5),
            # nn.Linear(256, 5)
            nn.Linear(256, 1)
        )

        self.base_channel = nn.Sequential(
            nn.Conv2d(1, 64, (6, 4)),
            nn.ReLU(inplace=True),
            nn.Flatten(start_dim=2),
            nn.MaxPool1d(3),
            nn.BatchNorm1d(64),
            nn.Conv1d(64, 128, 3),
            nn.ReLU(inplace=True),
            nn.MaxPool1d(3),
            nn.BatchNorm1d(128),
            nn.Conv1d(128, 256, 3),
            nn.ReLU(inplace=True),
            nn.AdaptiveAvgPool1d(1),
            nn.Flatten(),
        )

        self.codon_channel = nn.Sequential(
            nn.Conv2d(1, 64, (6, 64)),
            nn.ReLU(inplace=True),
            nn.Flatten(start_dim=2),
            nn.MaxPool1d(3),
            nn.BatchNorm1d(64),
            nn.Conv1d(64, 128, 3),
            nn.ReLU(inplace=True),
            nn.MaxPool1d(3),
            nn.BatchNorm1d(128),
            nn.Conv1d(128, 256, 3),
            nn.ReLU(inplace=True),
            nn.AdaptiveAvgPool1d(1),
            nn.Flatten(),
        )

        self.dropout = nn.Dropout(p=0.5, inplace=True)

    def forward(self, x, tfidf):
        interm_forward = self.base_channel(x)
        codon = self.codon_transformer(x)
        codon = self.codon_channel(codon)
        rev = torch.flip(x, dims=[-1, -2])
        interm_backward = self.base_channel(rev)
        codon_backward = self.codon_transformer(rev)
        codon_backward = self.codon_channel(codon_backward)
        z = torch.cat((interm_forward, codon, interm_backward, codon_backward, tfidf), dim=1)
        z = self.dropout(z)
        z = self.fc(z)
        return z