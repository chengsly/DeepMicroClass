from Bio import SeqIO
import os

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--input', dest='input')
parser.add_argument('--output_dir', dest='output_dir')

args = parser.parse_args()

# data_dir = os.path.join('..', 'data')

# source_fasta = os.path.join(data_dir, 'virushostdb.formatted.genomic.fna')
source_fasta = args.input

"""
Splitting virus genome into separated fasta files
"""
# split_dir = os.path.join(args.output_dir, 'virus_fasta')
split_dir = args.output_dir
os.makedirs(split_dir, exist_ok=True)
for record in SeqIO.parse(source_fasta, 'fasta'):
    # refseqid = record.description.split(' ')[0]
    refseqid = record.id.replace('/','')
    if not os.path.exists(os.path.join(split_dir, f'{refseqid}.fasta')):
        SeqIO.write(record, os.path.join(split_dir, f'{refseqid}.fasta'), 'fasta')
    