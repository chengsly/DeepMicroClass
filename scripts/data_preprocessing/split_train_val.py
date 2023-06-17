import argparse
import math
from Bio import SeqIO
import pandas as pd
import os
from ete3 import NCBITaxa
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("--input", "-i", dest="input")
parser.add_argument("--output_prefix", "-o", dest="output_prefix")
parser.add_argument("--dist", "-d", dest="dist")
parser.add_argument("--train_ratio", dest="train_ratio", default=0.8, type=float)
parser.add_argument("--threshold", dest="threshold", default=0.01, type=float)
args = parser.parse_args()

dist = pd.read_table(args.dist, index_col=0, delimiter="\t", engine="c", na_filter=False, low_memory=False)

org_records = list(SeqIO.parse(args.input, "fasta"))
org_id = [record.id for record in org_records]

record_number = len(org_id)

shuffled_index = np.arange(record_number)
np.random.shuffle(shuffled_index)

boundary = math.floor(record_number * args.train_ratio)

train_id = [org_id[i] for i in shuffled_index[:boundary]]
val_id = [org_id[i] for i in shuffled_index[boundary:]]

train_val_dist = dist.loc[train_id, val_id]
train_val_overlap = set(train_val_dist.columns[(train_val_dist < args.threshold).sum(axis=0) == 0])

train_id.extend(train_val_overlap)
train_id = set(train_id)
val_id = set(id for id in val_id if id not in train_val_overlap)

train_records = (record for record in org_records if record.id in train_id)
val_records = (record for record in org_records if record.id in val_id)

train_path = f"{args.output_prefix}Tr.fasta"
SeqIO.write(train_records, train_path, "fasta")

val_path = f"{args.output_prefix}Val.fasta"
SeqIO.write(val_records, val_path, "fasta")
