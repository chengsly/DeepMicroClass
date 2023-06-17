import argparse
import math
from Bio import SeqIO
import pandas as pd
import os
from ete3 import NCBITaxa
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("--train", "-i1", dest="train")
parser.add_argument("--test", "-i2", dest="test")
parser.add_argument("--output", "-o", dest="output")
parser.add_argument("--dist", "-d", dest="dist")
parser.add_argument("--threshold", dest="threshold", default=0.01, type=float)
args = parser.parse_args()

dist = pd.read_table(args.dist, index_col=0, delimiter='\t', engine='c', na_filter=False, low_memory=False)

train_id = [record.id for record in SeqIO.parse(args.train, 'fasta') if record.id in dist.columns]

test_records = list(SeqIO.parse(args.test, 'fasta'))
test_id = [record.id for record in test_records if record.id in dist.index]

# train_test_dist = dist.loc[train_id, test_id]
train_test_dist = dist.loc[test_id, train_id].T
train_test_overlap = set(train_test_dist.columns[(train_test_dist < args.threshold).sum(axis=0) != 0]) | set(train_test_dist.index[(train_test_dist < args.threshold).sum(axis=1) != 0])

test_id = set(id for id in test_id if id not in train_test_overlap)

test_records = (record for record in test_records if record.id in test_id)

SeqIO.write(test_records, args.output, 'fasta')