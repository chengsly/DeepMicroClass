from Bio import SeqIO
import pandas as pd
import os
from ete3 import NCBITaxa

DB_METADATA_PATH = "data/IMGVR/IMGVR_all_Sequence_information.tsv"
PROK_VIR_DIR = "data/IMGVR/prok_vir"
EUK_VIR_DIR = "data/IMGVR/euk_vir"
FASTA_PATH = "data/IMGVR/IMGVR_all_nucleotides.fna"

metadata = pd.read_table(DB_METADATA_PATH)
metadata.dropna(subset=["Host taxonomy prediction"], inplace=True)
metadata.reset_index(drop=True, inplace=True)

prok_vir = set()
euk_vir = set()
for i, row in metadata.iterrows():
    if row["Host taxonomy prediction"].split(";")[0] in {"Archaea", "Bacteria"}:
        prok_vir.add(row["## UViG"])
    elif row["Host taxonomy prediction"].split(";")[0] == "Eukaryota":
        euk_vir.add(row["## UViG"])

os.makedirs(PROK_VIR_DIR, exist_ok=True)
os.makedirs(EUK_VIR_DIR, exist_ok=True)

for record in SeqIO.parse(FASTA_PATH, "fasta"):
    uvig = record.id.split("|")[0]
    if uvig in prok_vir:
        directory = PROK_VIR_DIR
    elif uvig in euk_vir:
        directory = EUK_VIR_DIR
    else:
        continue
    output_path = os.path.join(directory, f"{record.id}.fa")
    SeqIO.write(record, output_path, "fasta")
