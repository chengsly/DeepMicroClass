from Bio import SeqIO
import pandas as pd
import os
from ete3 import NCBITaxa

DB_METADATA_PATH = "data/vhdb/virushostdb.tsv"
PROK_VIR_DIR = "data/vhdb/prok_vir"
EUK_VIR_DIR = "data/vhdb/euk_vir"
VHDB_FASTA_PATH = "data/vhdb/virushostdb.formatted.genomic.fna"

metadata = pd.read_table(DB_METADATA_PATH)
metadata.dropna(subset=["host tax id", "host lineage"], inplace=True)
metadata.reset_index(drop=True, inplace=True)

prok_vir = set()
euk_vir = set()
for i, row in metadata.iterrows():
    if row["host lineage"].split(";")[0] in {"Archaea", "Bacteria"}:
        for virus in row["refseq id"].split(", "):
            prok_vir.add(virus)
    elif row["host lineage"].split(";")[0] == "Eukaryota":
        for virus in row["refseq id"].split(", "):
            euk_vir.add(virus)

os.makedirs(PROK_VIR_DIR, exist_ok=True)
os.makedirs(EUK_VIR_DIR, exist_ok=True)

for record in SeqIO.parse(VHDB_FASTA_PATH, "fasta"):
    if record.id in prok_vir:
        directory = PROK_VIR_DIR
    elif record.id in euk_vir:
        directory = EUK_VIR_DIR
    output_path = os.path.join(directory, f"{record.id}.fa")
    SeqIO.write(record, output_path, "fasta")
