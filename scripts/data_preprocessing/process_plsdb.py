from Bio import SeqIO
import pandas as pd
import os
from ete3 import NCBITaxa
import numpy as np

DB_METADATA_PATH = "data/plsdb/plsdb.tsv"
FA_PATH = "data/plsdb/plsdb.fna"
PRE20_PATH = "data/plasmid/pre20_plasmid.fa"
POST20_PATH = "data/plasmid/post20_plasmid.fa"

sep_date = np.datetime64('2020-01-01')

metadata = pd.read_table(DB_METADATA_PATH)
dates = pd.to_datetime(metadata['SubmissionDate_ASSEMBLY'])
pre20_metadata = metadata.loc[dates < sep_date]
post20_metadata = metadata.loc[dates >= sep_date]

pre20_set = set(pre20_metadata['ACC_NUCCORE'])
post20_set = set(post20_metadata['ACC_NUCCORE'])

pre20_records = []
post20_records = []

for record in SeqIO.parse(FA_PATH, 'fasta'):
    if record.id in pre20_set:
        pre20_records.append(record)
    elif record.id in post20_set:
        post20_records.append(record)

SeqIO.write(pre20_records, PRE20_PATH, 'fasta')
SeqIO.write(post20_records, POST20_PATH, 'fasta')