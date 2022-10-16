import argparse
from Bio import SeqIO
import pandas as pd
import os
from ete3 import NCBITaxa
import numpy as np
from tqdm import tqdm

INPUT_DIR = 'data/sampled_sequence'
OUTPUT_DIR = 'data/sampled_sequence_subsampled_2000'
os.makedirs(OUTPUT_DIR, exist_ok=True)

parser = argparse.ArgumentParser()
parser.add_argument("--length", "-l", dest="length", default=10000, type=int, help="length of the sampled sequence")
args = parser.parse_args()

def sample_seq(seq, length):
    if len(seq) <= length:
        return seq
    start = np.random.randint(0, len(seq) - length)
    return seq[start:start+length]

fn = os.listdir(INPUT_DIR)
for f in tqdm(fn):
    input_generator = SeqIO.parse(os.path.join(INPUT_DIR, f), 'fasta')
    with open(os.path.join(OUTPUT_DIR, f), 'w') as output_file:
        output_record = [sample_seq(seq, int(args.length)) for seq in input_generator]
        SeqIO.write(output_record, output_file, 'fasta')