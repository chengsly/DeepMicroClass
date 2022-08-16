from Bio import SeqIO
import argparse
import gzip
import os
from tqdm import tqdm

parser = argparse.ArgumentParser()
parser.add_argument("--input_dir", "-i", dest="input_dir")
parser.add_argument("--output_dir", "-o", dest="output_dir")
parser.add_argument("--gzip", "-z", dest="use_gzip")

args = parser.parse_args()

os.makedirs(args.output_dir, exist_ok=True)

input_list = [fn for fn in os.listdir(args.input_dir) if fn.endswith(("gz", "fasta", "fna", "fa"))]

for fasta_file in tqdm(input_list):
    input_fasta_path = os.path.join(args.input_dir, fasta_file)
    if args.use_gzip:
        input_handle = gzip.open(input_fasta_path, "rt")
    else:
        input_handle = open(input_fasta_path, "r")

    output_fasta_path = os.path.join(args.output_dir, fasta_file)
    if args.use_gzip:
        output_handle = gzip.open(output_fasta_path, "wt")
    else:
        output_handle = open(output_fasta_path, "w")

    seq_iter = SeqIO.parse(input_handle, format="fasta")
    seq_list = []
    for seq in seq_iter:
        if "plasmid" not in seq.description.lower():
            seq_list.append(seq)
    SeqIO.write(seq_list, output_handle, "fasta")

    input_handle.close()
    output_handle.close()
