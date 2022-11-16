from Bio import SeqIO
import pandas as pd
import os
from ete3 import NCBITaxa

DB_METADATA_PATH = "/home/tianqi/dataset/IMGVR_v2/IMGVR_all_Sequence_information.tsv"
OUTPUT_PATH = "/home/tianqi/dataset/IMGVR_v2/sources"
FASTA_PATH = "/home/tianqi/dataset/IMGVR_v2/IMGVR_all_nucleotides.fna"

# DB_METADATA_PATH = "/home/tianqi/dataset/IMGVR_v4_high_confidence/IMGVR_all_Sequence_information.tsv"
# OUTPUT_PATH = "/home/tianqi/dataset/IMGVR_v4_high_confidence/sources"
# FASTA_PATH = "/home/tianqi/dataset/IMGVR_v4_high_confidence/IMGVR_all_nucleotides.fna"

os.makedirs(OUTPUT_PATH, exist_ok=True)

metadata = pd.read_table(DB_METADATA_PATH)
metadata.dropna(subset=["Ecosystem_Category", "Ecosystem_Type", "Ecosystem_Subtype"], inplace=True)
# metadata.dropna(subset=["Ecosystem classification"], inplace=True)
metadata.reset_index(drop=True, inplace=True)

conds = ['Animal', 'Aquatic_Sediment', 'Plants']
# conds_sets = [set(metadata[metadata["Ecosystem classification"].str.contains(cond)]["## UViG"]) for cond in conds]
# conds_sets = [set(metadata[metadata["Ecosystem_Category"].str.contains(cond)]["UViG"]) for cond in conds]
# conds_sets = [set(metadata[metadata["Ecosystem classification"].str.contains(cond)]["UVIG"]) for cond in conds]
conds_sets = [
    set(metadata[metadata["Ecosystem_Category"].str.contains('Animal')]["UViG"]),
    set(metadata[metadata["Ecosystem_Type"].str.contains('Sediment')]["UViG"]),
    set(metadata[metadata["Ecosystem_Category"].str.contains('Plants')]["UViG"]),
]
handles = [open(os.path.join(OUTPUT_PATH, f"{cond}.fasta"), "w") for cond in conds]

for record in SeqIO.parse(FASTA_PATH, "fasta"):
    uvig = record.id.split("|")[0]
    for i, cond in enumerate(conds):
        if uvig in conds_sets[i]:
            SeqIO.write(record, handles[i], "fasta")
    # output_path = os.path.join(directory, f"{record.id}.fa")
    # SeqIO.write(record, output_path, "fasta")
