import numpy as np
import pandas as pd
import os
import re
from Bio import SeqIO

FASTA_DIR = 'data/sampled_sequence'

results_fn = os.listdir(FASTA_DIR)
results_fn.sort()

TARGET_DIR = 'data/sampled_sequence_target'
os.makedirs(TARGET_DIR, exist_ok=True)

def construct_target(nums):
    order = [[3], [0], [2], [4], [1]]
    return np.array(sum([order[i]*nums[i] for i in range(5)], []))

for f in results_fn:
    nums = re.findall(r'\d+', f)[1:]
    nums = [int(n) for n in nums]
    target = construct_target(nums)

    ids = [seq.id for seq in SeqIO.parse(os.path.join(FASTA_DIR, f), 'fasta')]

    df = pd.DataFrame(target, index=ids)

    df.to_pickle(os.path.join(TARGET_DIR, f'{f}.pkl'))