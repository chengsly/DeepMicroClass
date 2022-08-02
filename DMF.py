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

class DMF(nn.Module):
    def __init__(self) -> None:
        super().__init__()

        self.codon_transformer = CodonTransformer()

        self.fc = nn.Sequential(
            # nn.LazyLinear(512),
            # nn.ReLU(inplace=True),
            # nn.Dropout(p=0.1),
            # nn.Linear(512, 256),
            nn.Linear(1024, 256),
            nn.PReLU(),
            nn.Dropout(p=0.2),
            nn.Linear(256, 5)
        )

        self.base_channel = nn.Sequential(
            nn.Conv2d(1, 64, (6, 4)),
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
        z = self.dropout(z)
        z = self.fc(z)
        return z

class DMF_3class(nn.Module):
    def __init__(self) -> None:
        super().__init__()

        self.codon_transformer = CodonTransformer()

        self.fc = nn.Sequential(
            nn.Linear(512, 256),
            nn.PReLU(),
            nn.Dropout(p=0.2),
            nn.Linear(256, 3)
        )

        self.base_channel = nn.Sequential(
            # nn.Conv2d(1, 64, (6, 4)),
            # nn.PReLU(),
            # nn.Flatten(start_dim=2),
            # # nn.MaxPool1d(3),
            # nn.AvgPool1d(3),
            # nn.BatchNorm1d(64),
            # nn.Conv1d(64, 128, 3),
            # nn.PReLU(),
            # # nn.MaxPool1d(3),
            # nn.AvgPool1d(3),
            # nn.BatchNorm1d(128),
            # nn.Conv1d(128, 256, 3),
            # nn.PReLU(),
            # nn.BatchNorm1d(256),
            # nn.AdaptiveAvgPool1d(1),
            # nn.Flatten(),

            nn.Conv1d(4, 64, 6),
            nn.PReLU(),
            nn.MaxPool1d(stride=2, kernel_size=2),
            nn.BatchNorm1d(64),
            nn.Conv1d(64, 128, 3),
            nn.PReLU(),
            nn.MaxPool1d(stride=1, kernel_size=2),
            nn.BatchNorm1d(128),
            nn.Conv1d(128, 256, 2),
            nn.PReLU(),
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
            nn.AdaptiveAvgPool1d(1),
            nn.Flatten(),
        )

        self.dropout = nn.Dropout(p=0.2, inplace=True)

    def forward(self, x):
        x = x.permute(0, 2, 1)
        interm_forward = self.base_channel(x)
        # codon = self.codon_transformer(x)
        # codon = self.codon_channel(codon)
        rev = torch.flip(x, dims=[-1, -2])
        interm_backward = self.base_channel(rev)
        # codon_backward = self.codon_transformer(rev)
        # codon_backward = self.codon_channel(codon_backward)
        # z = torch.cat((interm_forward, codon, interm_backward, codon_backward, ), dim=1)
        z = torch.cat((interm_forward, interm_backward, ), dim=1)
        z = self.dropout(z)
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

class LightningDMF(pl.LightningModule):
    def __init__(self, model, weight=None, num_classes=5):
        super().__init__()
        self.model = model
        self.weight = weight
        self.f1_log_handle = open("log/f1_log.txt", "w")
        self.num_classes = num_classes

    def forward(self, x):
        return self.model(x)

    def training_step(self, batch, batch_idx):
        x, y = batch
        y_hat = self.model(x)
        loss = F.cross_entropy(y_hat, y)
        # loss = F.cross_entropy(y_hat, y, weight=self.weight)
        # loss = F.binary_cross_entropy_with_logits(y_hat, y, pos_weight=self.weight[1])
        self.log("train_loss", loss)
        # self.log('train_acc_step', accuracy(y_hat, y.int(), average='macro'), on_step=True, on_epoch=False, prog_bar=True)
        self.log('train_acc_epoch', accuracy(y_hat, y.int(), average='macro', num_classes=self.num_classes), on_step=False, on_epoch=True, prog_bar=True)
        self.log('train_f1_epoch', f1_score(y_hat, y.int(), average='macro', num_classes=self.num_classes), on_step=False, on_epoch=True, prog_bar=True)
        return loss
    
    # def training_epoch_end(self, outputs) -> None:
    #     avg_acc = torch.stack([x["train_acc_step"] for x in outputs]).mean()
    #     self.log("train_acc_epoch", avg_acc, prog_bar=True)
    #     return super().training_epoch_end(outputs)

    def validation_step(self, batch, batch_idx):
        x, y = batch
        y_hat = self.model(x)
        loss = F.cross_entropy(y_hat, y)
        # loss = F.cross_entropy(y_hat, y, weight=self.weight)
        # loss = F.binary_cross_entropy_with_logits(y_hat, y, pos_weight=self.weight[1])
        self.log("val_loss", loss, prog_bar=True)
        self.log('val_acc', accuracy(y_hat, y.int(), average='macro', num_classes=self.num_classes), prog_bar=True)
        self.log('val_f1', f1_score(y_hat, y.int(), average='macro', num_classes=self.num_classes), prog_bar=True)
        # print(f1_score(y_hat, y.int(), average=None, num_classes=5), end='\r')
        self.f1_log_handle.write(f'{f1_score(y_hat, y.int(), average=None, num_classes=self.num_classes)}\n')
        # self.log('val_f1_plasmid', f1_score(y_hat, y.int(), average=None, num_classes=5)[2], prog_bar=False, on_step=False, on_epoch=True)
        return loss
    
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
        return torch.optim.Adam(self.parameters(), lr=1e-3)

class DMF_2class(nn.Module):
    def __init__(self) -> None:
        super().__init__()

        self.codon_transformer = CodonTransformer()

        self.fc = nn.Sequential(
            nn.LazyLinear(512),
            nn.ReLU(inplace=True),
            nn.Dropout(p=0.1),
            nn.Linear(512, 256),
            nn.ReLU(inplace=True),
            nn.Dropout(p=0.1),
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

        self.dropout = nn.Dropout(p=0.1, inplace=True)

    def forward(self, x):
        interm_forward = self.base_channel(x)
        codon = self.codon_transformer(x)
        codon = self.codon_channel(codon)
        rev = torch.flip(x, dims=[-1, -2])
        interm_backward = self.base_channel(rev)
        codon_backward = self.codon_transformer(rev)
        codon_backward = self.codon_channel(codon_backward)
        z = torch.cat((interm_forward, codon, interm_backward, codon_backward, ), dim=1)
        z = self.dropout(z)
        z = self.fc(z)
        return z

class LightningDMF_2class(pl.LightningModule):
    def __init__(self, model, weight=None):
        super().__init__()
        self.model = model
        self.weight = weight

    def forward(self, x):
        return self.model(x)

    def training_step(self, batch, batch_idx):
        x, y = batch
        y_hat = self.model(x).flatten()
        # loss = F.cross_entropy(y_hat, y)
        # loss = F.cross_entropy(y_hat, y, weight=self.weight)
        loss = F.binary_cross_entropy_with_logits(y_hat, y)
        self.log("train_loss", loss)
        # self.log('train_acc_epoch', accuracy(F.sigmoid(y_hat), y.int(), num_classes=1), on_step=False, on_epoch=True, prog_bar=True)
        self.log('train_f1_epoch', f1_score(F.sigmoid(y_hat), y.int()), on_step=False, on_epoch=True, prog_bar=True)
        return loss

    def validation_step(self, batch, batch_idx):
        x, y = batch
        y_hat = self.model(x).flatten()
        loss = F.binary_cross_entropy_with_logits(y_hat, y)
        self.log("val_loss", loss, prog_bar=True)
        # self.log('val_acc', accuracy(F.sigmoid(y_hat), y.int(), num_classes=1), prog_bar=True)
        self.log('val_f1', f1_score(F.sigmoid(y_hat), y.int()), prog_bar=True)
        return loss

    def configure_optimizers(self):
        return torch.optim.Adam(self.parameters(), lr=1e-4)

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

class DMFTransformer(nn.Module):
    def __init__(self) -> None:
        super().__init__()

        self.kmer_transformer = KMerTransformer(k=4)

        # self.embeddings = [nn.Embedding(4 ** 4, 32), nn.Embedding(4 ** 5, 128), nn.Embedding(4 ** 6, 512)]
        # self.pos_encoders = [PositionalEncoding(32, 0.1), PositionalEncoding(128, 0.1), PositionalEncoding(512, 0.1)]
        # self.embedding_5 = nn.Embedding(4 ** 5, 128)
        # self.embedding_6 = nn.Embedding(4 ** 6, 512)
        self.k = 5
        self.embed_d = 128
        self.embedding = nn.Embedding(4 ** 4, self.embed_d)
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
    
    def forward(self, src):
        src = self.kmer_transformer(src)
        src = torch.argmax(src, dim=1).long()
        src = self.embedding(src) * math.sqrt(self.embed_d)
        src = self.pos_encoder(src)
        output = self.transformer_encoder(src)
        output = output.mean(dim=1)
        output = self.decoder(output)
        return output

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