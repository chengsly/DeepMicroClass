from Bio import SeqIO
from Bio import Entrez
import pandas as pd
import os
from ete3 import NCBITaxa

DB_METADATA_PATH = "data/vhdb/virushostdb.tsv"
PROK_VIR_DIR = "data/vhdb/prok_vir"
EUK_VIR_DIR = "data/vhdb/euk_vir"
VHDB_FASTA_PATH = "data/vhdb/virushostdb.formatted.genomic.fna"
UPDATE_DATES_PATH = "data/vhdb/update_dates.tsv"
PROK_PRE20_PATH = "data/vhdb/prokvir_pre20.fa"
PROK_POST20_PATH = "data/vhdb/prokvir_post20.fa"
EUK_PRE20_PATH = "data/vhdb/eukvir_pre20.fa"
EUK_POST20_PATH = "data/vhdb/eukvir_post20.fa"

if os.path.exists(UPDATE_DATES_PATH):
    updates_dates = pd.read_table(UPDATE_DATES_PATH, sep="\t", header=None, names=["refseq id", "update date"], index_col=0)
    # to dict
    updates_dates = updates_dates.to_dict()["update date"]
else:
    updates_dates = {}

update_date_handle = open(UPDATE_DATES_PATH, "w")
for refseq_id, update_date in updates_dates.items():
    update_date_handle.write(f"{refseq_id}\t{update_date}\n")
update_date_handle.flush()

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

# os.makedirs(PROK_VIR_DIR, exist_ok=True)
# os.makedirs(EUK_VIR_DIR, exist_ok=True)

# for record in SeqIO.parse(VHDB_FASTA_PATH, "fasta"):
#     if record.id in prok_vir:
#         directory = PROK_VIR_DIR
#     elif record.id in euk_vir:
#         directory = EUK_VIR_DIR
#     output_path = os.path.join(directory, f"{record.id}.fa")
#     SeqIO.write(record, output_path, "fasta")

# Separate viruses according to 2020/01/01 update date
def get_update_date(refseq_id):
    Entrez.email = "tianqit@usc.edu"
    handle = Entrez.esummary(db="nucleotide", id=refseq_id)
    record = Entrez.read(handle)
    handle.close()
    print(f"{refseq_id}: {record[0]['UpdateDate']}", flush=True, end="\r")
    return record[0]["UpdateDate"]

prok_pre20 = []
prok_post20 = []
euk_pre20 = []
euk_post20 = []

for record in SeqIO.parse(VHDB_FASTA_PATH, "fasta"):
    if record.id in prok_vir:
        if record.id in updates_dates:
            update_date = updates_dates[record.id]
        else:
            update_date = get_update_date(record.id)
            updates_dates[record.id] = update_date
            # save updates_dates to tsv
            update_date_handle.write(f"{record.id}\t{update_date}\n")
            update_date_handle.flush()
        if update_date < "2020/01/01":
            prok_pre20.append(record)
        else:
            prok_post20.append(record)
    elif record.id in euk_vir:
        if record.id in updates_dates:
            update_date = updates_dates[record.id]
        else:
            update_date = get_update_date(record.id)
            updates_dates[record.id] = update_date
            # save updates_dates to tsv
            update_date_handle.write(f"{record.id}\t{update_date}\n")
            update_date_handle.flush()
        if update_date < "2020/01/01":
            euk_pre20.append(record)
        else:
            euk_post20.append(record)

SeqIO.write(prok_pre20, PROK_PRE20_PATH, "fasta")
SeqIO.write(prok_post20, PROK_POST20_PATH, "fasta")
SeqIO.write(euk_pre20, EUK_PRE20_PATH, "fasta")
SeqIO.write(euk_post20, EUK_POST20_PATH, "fasta")
