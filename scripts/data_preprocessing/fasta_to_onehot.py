import sys
sys.path.append('..')

import numpy as np
import util
import os
from Bio import SeqIO, Data
from tqdm import tqdm

from DeepMicroClass.constants import CODON_TABLE

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--input_dir', dest='input_dir')
parser.add_argument('--output_dir', dest='output_dir')
parser.add_argument('--plasmid', dest='plasmid', action='store_true')

args = parser.parse_args()

os.makedirs(args.output_dir, exist_ok=True)
os.makedirs(args.output_dir+'_aa', exist_ok=True)

fasta_ext = ('.fna', '.fasta', '.fa')

input_fasta = os.listdir(args.input_dir)
input_fasta = list(filter(lambda x: x.endswith(fasta_ext), input_fasta))

possible_aa = set(Data.IUPACData.protein_letters)
possible_aa.add('*')
possible_aa.add('X')
possible_aa.add('B')
possible_aa.add('Z')
possible_aa.add('J')
possible_aa.add('O')
possible_aa.add('U')
aa_index = {aa: i for i, aa in enumerate(possible_aa)}

for fasta in tqdm(input_fasta):
    fasta_path = os.path.join(args.input_dir, fasta)

    output_onehot_fn = fasta + '.npy'
    output_onehot_path = os.path.join(args.output_dir, output_onehot_fn)

    # if os.path.exists(output_onehot_path):
    #     continue
    
    seq_records = list(SeqIO.parse(fasta_path, 'fasta'))
    if not os.path.exists(output_onehot_path):
        if not args.plasmid:
            seq = [str(record.seq) for record in seq_records if not 'plasmid' in record.description]
        else:
            seq = [str(record.seq) for record in seq_records]
        seq = ''.join(seq)

        onehot = util.seq2onehot(seq)
        onehot = onehot[:, :4]
        
        np.save(output_onehot_path, onehot.astype(np.int8))

    aa_sequences = [record.translate() for record in seq_records]
    aa_sequences1 = [record[1:].translate() for record in seq_records]
    aa_sequences2 = [record[2:].translate() for record in seq_records]
    if not args.plasmid:
        aa_sequences = [str(record.seq) for record in aa_sequences if not 'plasmid' in record.description]
        aa_sequences1 = [str(record.seq) for record in aa_sequences1 if not 'plasmid' in record.description]
        aa_sequences2 = [str(record.seq) for record in aa_sequences2 if not 'plasmid' in record.description]
    else:
        aa_sequences = [str(record.seq) for record in aa_sequences]
        aa_sequences1 = [str(record.seq) for record in aa_sequences1]
        aa_sequences2 = [str(record.seq) for record in aa_sequences2]
    aa_sequences = ''.join(aa_sequences + aa_sequences1 + aa_sequences2)

    aa_onehot = np.zeros((len(aa_sequences), 27), dtype=np.int8)
    aa_indices = [aa_index[aa] for aa in aa_sequences]
    aa_onehot[np.arange(len(aa_indices)), aa_indices] = 1

    output_aa_path = os.path.join(args.output_dir+'_aa', output_onehot_fn)
    np.save(output_onehot_path + '.aa.npy', aa_onehot.astype(np.int8))