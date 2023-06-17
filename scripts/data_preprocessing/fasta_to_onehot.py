import sys
sys.path.append('..')

import numpy as np
import util
import os
from Bio import SeqIO
from tqdm import tqdm

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--input_dir', dest='input_dir')
parser.add_argument('--output_dir', dest='output_dir')

args = parser.parse_args()

os.makedirs(args.output_dir, exist_ok=True)

fasta_ext = ('.fna', '.fasta', '.fa')

input_fasta = os.listdir(args.input_dir)
input_fasta = list(filter(lambda x: x.endswith(fasta_ext), input_fasta))

for fasta in tqdm(input_fasta):
    fasta_path = os.path.join(args.input_dir, fasta)

    output_onehot_fn = fasta + '.npy'
    output_onehot_path = os.path.join(args.output_dir, output_onehot_fn)

    if os.path.exists(output_onehot_path):
        continue

    seq_records = list(SeqIO.parse(fasta_path, 'fasta'))
    # seq = [str(record.seq) for record in seq_records if not 'plasmid' in record.description]
    seq = [str(record.seq) for record in seq_records]
    seq = ''.join(seq)

    onehot = util.seq2onehot(seq)
    onehot = onehot[:, :4]

    
    
    np.save(output_onehot_path, onehot.astype(np.int8))