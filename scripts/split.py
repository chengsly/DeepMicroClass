import argparse
import math
from Bio import SeqIO
import pandas as pd
import os
from ete3 import NCBITaxa
import numpy as np
from random import shuffle

parser = argparse.ArgumentParser()
parser.add_argument("--input", "-i", dest="input")
parser.add_argument("--output_prefix", "-o", dest="output_prefix")
parser.add_argument("--train_ratio", dest="train_ratio", default=0.7, type=float)
parser.add_argument("--max-seq-n", "-m", dest="max_seq", type=int)
args = parser.parse_args()

org_records = list(SeqIO.parse(args.input, "fasta"))
org_id = [record.id for record in org_records]

shuffle(org_id)
org_id = org_id[:args.max_seq]

border = math.floor(len(org_id) * args.train_ratio)

train_id = set(org_id[:border])
val_id = set(org_id[border:])

train_records = (record for record in org_records if record.id in train_id)
val_records = (record for record in org_records if record.id in val_id)

train_path = f"{args.output_prefix}Tr.fasta"
SeqIO.write(train_records, train_path, "fasta")

val_path = f"{args.output_prefix}Val.fasta"
SeqIO.write(val_records, val_path, "fasta")